import numpy as np
import pandas as pd
import os
import pickle
from functools import reduce
import anndata as ad
import multiprocessing as mp
import pysam
import warnings
from pathlib import Path
from sklearn.cluster import AgglomerativeClustering
from scipy.optimize import linear_sum_assignment
import time
import subprocess
import re
import sys
from .deep_learning import scDataset, Net,test 
from scipy import stats
from torch.utils.data import DataLoader
from .toolbox import check_pysam_chrom,fetch_reads
from sklearn.neighbors import KernelDensity
import racplusplus
from collections import Counter
import matplotlib.pyplot as plt
import random 
from scipy.sparse import csr_matrix
from .build_ref import get_gene_with_one_transcript


warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=Warning)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



class GetPAScluster():
    def __init__(self,PASrefPath,generefPath,fastafilePath,bamfilePath,outdir,nproc,minCount,maxReadCount,densityFC,InnerDistance,device,chromoSizePath,cellbarcodePath,species):
    
        self.PASrefdf=pd.read_csv(PASrefPath,delimiter='\t')
        self.generefdf=pd.read_csv(generefPath,delimiter='\t')
        self.chromosizedf=pd.read_csv(chromoSizePath,delimiter='\t',header=None,index_col=0)
        self.cellbarcode=pd.read_csv(cellbarcodePath,delimiter='\t')['cellbarcode'].values
        self.fastafilePath=fastafilePath
        self.bamfilePath=bamfilePath
        self.outdir=outdir
        self.count_out_dir=os.path.join(outdir,'count')
        if not os.path.exists(self.count_out_dir):
            os.mkdir(self.count_out_dir)

        self.filteredbamfilePath=os.path.join(self.count_out_dir,'filtered_reads.bam')
        self.pcr_removedPath=os.path.join(self.count_out_dir,'pcr_duplication_removed.bam')

        self.minCount=minCount
        self.nproc=nproc
        self.maxReadCount=maxReadCount
        self.densityFC=densityFC
        self.InnerDistance=InnerDistance
        self.device=device
        self.species=species 




    def _getreads(self,bamFile,chr):


        samFile, _chrom = check_pysam_chrom(bamFile, str(chr))
        reads = fetch_reads(samFile, _chrom,  0 , self.chromosizedf.loc[chr][1],  trimLen_max=100)
        reads1 = reads["reads1"]   

        reads1=[r for r in reads1 if r.get_tag('CB') in self.cellbarcode]


        forwardReads1=[r for r in reads1 if r.is_reverse==True]
        reverseReads1=[r for r in reads1 if r.is_reverse==False]




        #customized forward remove pcr duplication; want to use the reference_end whose value is the largest. 10x genomics, first do pcr amplification and then do fragmentation. 
        #So just require the cellbarcode, UMI are the same, even the sequence is not the same, the alignment site is not the same, they still are the duplication. 
        forwardreadsdf=pd.DataFrame({'CB':[r.get_tag('CB') for r in forwardReads1],'UMI':[r.get_tag('UB') for r in forwardReads1],'gene':[r.get_tag('GX') for r in forwardReads1],'PAS':[r.reference_end for r in forwardReads1]})
        notoperate=forwardreadsdf[~forwardreadsdf.duplicated(['CB','UMI','gene'],keep=False)]
        duplicationdf=forwardreadsdf[forwardreadsdf.duplicated(['CB','UMI','gene'],keep=False)]
        duplicationdf.sort_values(['CB','UMI','gene','PAS'],inplace=True)
        afterduplication=duplicationdf.groupby(['CB','UMI','gene']).tail(1)
        forward_pcr_deduplicationdf=pd.concat([notoperate,afterduplication],axis=0)
        reads1_PAS_forward=forward_pcr_deduplicationdf['PAS'].tolist()


        reversereadsdf=pd.DataFrame({'CB':[r.get_tag('CB') for r in reverseReads1],'UMI':[r.get_tag('UB') for r in reverseReads1],'gene':[r.get_tag('GX') for r in reverseReads1],'PAS':[r.reference_start for r in reverseReads1]})
        reverse_notoperate=reversereadsdf[~reversereadsdf.duplicated(['CB','UMI','gene'],keep=False)]
        reverse_duplicationdf=reversereadsdf[reversereadsdf.duplicated(['CB','UMI','gene'],keep=False)]
        reverse_duplicationdf.sort_values(['CB','UMI','gene','PAS'],inplace=True)
        reverse_afterduplication=reverse_duplicationdf.groupby(['CB','UMI','gene']).head(1)
        reverse_pcr_deduplicationdf=pd.concat([reverse_notoperate,reverse_afterduplication],axis=0)
        reads1_PAS_reverse=reverse_pcr_deduplicationdf['PAS'].tolist()



        forward_pos,forward_count=np.unique(reads1_PAS_forward,return_counts=True)
        reverse_pos,reverse_count=np.unique(reads1_PAS_reverse,return_counts=True)

        forwarddf=pd.DataFrame({'pos':forward_pos,'count':forward_count})
        forwarddf['strand']='+'

        reversedf=pd.DataFrame({'pos':reverse_pos,'count':reverse_count})
        reversedf['strand']='-'

        onechrdf=pd.concat([forwarddf,reversedf],axis=0)
        onechrdf['chr']=chr

        return onechrdf




    def _get_reads(self):

        start_time=time.time()
        bamFile=self.bamfilePath

        pool = mp.Pool(processes=self.nproc)


        allchrls=[]
        for chr in self.chromosizedf.index:
            allchrls.append(pool.apply_async(self._getreads,(bamFile,chr)))

        pool.close()
        pool.join()

        allchrls=[res.get() for res in allchrls]


        paracluinputdf=pd.concat(allchrls,axis=0)
        paracluinputdf=paracluinputdf[['chr','strand','pos','count']]
        paracluinputdf.sort_values(['chr','strand','pos'],inplace=True)
        paracluinputdf['pos']=paracluinputdf['pos'].astype('int')

        paraclu_inputPath=os.path.join(self.count_out_dir,'paraclu_input.tsv')
        paracluinputdf.to_csv(paraclu_inputPath,sep='\t',header=None,index=None)


        return paraclu_inputPath 




    
    def _do_cluster(self):
        start_time=time.time()

        paraclu_inputPath=self._get_reads()

        paraclu_outputPath=os.path.join(self.count_out_dir,'paraclu_output.tsv')

        #paraclu_CMD="paraclu {} {} > {}".format(self.minCount,paraclu_inputPath,paraclu_outputPath)
        paraclu_CMD="paraclu {} {} | paraclu-cut -l {} -d {} > {}".format(self.minCount, paraclu_inputPath, self.InnerDistance, self.densityFC, paraclu_outputPath)
        eprint(paraclu_CMD)
        subprocess.run(paraclu_CMD, shell=True,stdout=subprocess.PIPE)
        print('do clustering Time elapsed %.2f min' %((time.time()-start_time)/60))


        return paraclu_outputPath









    def _filter_false_positive(self):
        start_time=time.time()
        paraclu_outputPath=self._do_cluster()

        paracludf=pd.read_csv(paraclu_outputPath,header=None,delimiter='\t')
        paracludf['cluster_id']=paracludf[0]+'_'+paracludf[1]+'_'+paracludf[2].astype('str')+'_'+paracludf[3].astype('str')
        paracludf['PAS']=np.where(paracludf[1]=='+',paracludf[3],paracludf[2])
        inputtodeepdf=paracludf[[0,1,'PAS','cluster_id']]
        inputtodeepdf.columns=['Chromosome','Strand','PAS','cluster_id']



        input_to_DP=os.path.join(self.count_out_dir,'input_to_DP.tsv')
        inputtodeepdf.to_csv(input_to_DP,sep='\t',index=None)

        #run pre-trained deep learning model 
        test_set=scDataset(self.fastafilePath,input_to_DP,self.chromosizedf)
        test_loader = DataLoader(test_set, batch_size=1, shuffle=False)

        net=Net().to(self.device)

        if self.species=='human':
            pretrained_model_path=str(Path(os.path.dirname(os.path.abspath(__file__))).parents[0])+'/model/human_pretrained_model.pth'
        elif self.species=='mouse':
            pretrained_model_path=str(Path(os.path.dirname(os.path.abspath(__file__))).parents[0])+'/model/mouse_pretrained_model.pth'


        
        predict_result_output=os.path.join(self.count_out_dir,'predict_result.tsv')
        positivedf=test(test_loader,self.device,net,pretrained_model_path,predict_result_output)

        positive_result_output_path=os.path.join(self.count_out_dir,'positive_result.bed')

        positivedf['chr']=positivedf['PAS'].str.split('_',expand=True)[0]
        positivedf['strand']=positivedf['PAS'].str.split('_',expand=True)[1]
        positivedf['start']=positivedf['PAS'].str.split('_',expand=True)[2].astype('int')
        positivedf['end']=positivedf['PAS'].str.split('_',expand=True)[3].astype('int')
        positivedf['score']=0
        positivebeddf=positivedf[['chr','start','end','PAS','score','strand']]

        positivebeddf.to_csv(positive_result_output_path,sep='\t',header=None,index=None)

        print('Filtering false positive elapsed %.2f min' %((time.time()-start_time)/60))

        return positive_result_output_path







