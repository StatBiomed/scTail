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


warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=Warning)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



class get_PAS_count():
    def __init__(self,PASrefPath,generefPath,fastafilePath,bamfilePath,outdir,nproc,minCount,maxReadCount,clusterDistance,InnerDistance,device,chromoSizePath):
    
        self.PASrefdf=pd.read_csv(PASrefPath,delimiter='\t')
        self.generefdf=pd.read_csv(generefPath,delimiter='\t')
        self.chromosizedf=pd.read_csv(chromoSizePath,delimiter='\t',header=None,index_col=0)
        self.fastafilePath=fastafilePath
        self.bamfilePath=bamfilePath
        self.outdir=outdir
        self.count_out_dir=os.path.join(outdir,'count')
        if not os.path.exists(self.count_out_dir):
            os.mkdir(self.count_out_dir)

        self.filteredbamfilePath=os.path.join(self.count_out_dir,'filtered_reads.bam')
        self.pcr_removedPath=os.path.join(self.count_out_dir,'pcr_duplication_removed.bam')

        self.minCount=minCount
        #self.cellBarcode=pd.read_csv(cellBarcodePath,delimiter='\t')['cell_id'].values
        self.nproc=nproc
        self.maxReadCount=maxReadCount
        self.clusterDistance=clusterDistance
        self.InnerDistance=InnerDistance
        self.device=device





    # def _do_preprocess(self):
    #     start_time=time.time()

    #     filteredbamfilePath=self.filteredbamfilePath
    #     infile = pysam.AlignmentFile(self.bamfilePath, "rb")
    #     outfile = pysam.AlignmentFile(filteredbamfilePath, "wb", template=infile)
    #     for read in infile:
    #         if read.get_tag('GX')!='-':
    #             outfile.write(read)
    #     infile.close()
    #     outfile.close()

    #     pysam.index(self.filteredbamfilePath)
    #     print('filtering reads elapsed %.2f min' %((time.time()-start_time)/60))

    #     return filteredbamfilePath




    def _get_reads(self):

        start_time=time.time()
        bamFile=self.bamfilePath

        allchrls=[]
        for chr in self.chromosizedf.index:
            samFile, _chrom = check_pysam_chrom(bamFile, str(chr))
            reads = fetch_reads(samFile, _chrom,  0 , self.chromosizedf.loc[chr][1],  trimLen_max=100)
            reads1_umi = reads["reads2"]

            forwardReads1=[r for r in reads1_umi if r.is_reverse==True]
            reverseReads1=[r for r in reads1_umi if r.is_reverse==False]

            reads1_PAS_forward=[r.reference_end for r in forwardReads1]
            reads1_PAS_reverse=[r.reference_start for r in reverseReads1]

            forward_pos,forward_count=np.unique(reads1_PAS_forward,return_counts=True)
            reverse_pos,reverse_count=np.unique(reads1_PAS_reverse,return_counts=True)

            forwarddf=pd.DataFrame({'pos':forward_pos,'count':forward_count})
            forwarddf['strand']='+'

            reversedf=pd.DataFrame({'pos':reverse_pos,'count':reverse_count})
            reversedf['strand']='-'

            onechrdf=pd.concat([forwarddf,reversedf],axis=0)
            onechrdf['chr']=chr
            allchrls.append(onechrdf)


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
        paraclu_CMD="paraclu {} {} | paraclu-cut -l 100 -d 0 -s > {}".format(self.minCount, paraclu_inputPath, paraclu_outputPath)
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
        pretrained_model_path=str(Path(os.path.dirname(os.path.abspath(__file__))).parents[0])+'/model/model_last.pth'
        predict_result_output=os.path.join(self.count_out_dir,'predict_result.tsv')
        positivedf=test(test_loader,self.device,net,pretrained_model_path,predict_result_output)


        positivedf['Chromosome']=positivedf['PAS'].str.split('_',expand=True)[0]
        positivedf['Strand']=positivedf['PAS'].str.split('_',expand=True)[1]
        positivedf['Start']=positivedf['PAS'].str.split('_',expand=True)[2]
        positivedf['End']=positivedf['PAS'].str.split('_',expand=True)[3]
        featurecountinputdf=positivedf[['PAS','Chromosome','Start','End','Strand']]
        featurecountinputdf.columns=['GeneID','Chr','Start','End','Strand']

        featureCount_inputPath=os.path.join(self.count_out_dir,'featureCounts_input.tsv')
        featurecountinputdf.to_csv(featureCount_inputPath,sep='\t',index=None) 

        print('Filtering false positive elapsed %.2f min' %((time.time()-start_time)/60))

        return featureCount_inputPath



    def get_count_h5ad(self):

        ctime=time.time()

        featureCount_inputPath=self._filter_false_positive()
        featureCount_outputPath=os.path.join(self.count_out_dir,'featureCount_output')

        featureCount_CMD="featureCounts -a {} -o {} -F SAF -R BAM {} -T {}".format(featureCount_inputPath,featureCount_outputPath,self.pcr_removedPath,self.nproc)
        eprint(featureCount_CMD)
        subprocess.run(featureCount_CMD, shell=True,stdout=subprocess.PIPE)

        umi_tools_input=os.path.join(self.count_out_dir,'pcr_duplication_removed.bam.featureCounts.bam')
        sort_output=os.path.join(self.count_out_dir,'featurecount_sorted.bam')
        umi_tools_output=os.path.join(self.count_out_dir,'umitools_output')


        samtools_index_CMD="samtools sort -@ {} {} > {}".format(self.nproc,umi_tools_input,sort_output)
        eprint(samtools_index_CMD)
        subprocess.run(samtools_index_CMD,shell=True,stdout=subprocess.PIPE)


        samtools_index_CMD="samtools index -@ {} {}".format(self.nproc,sort_output)
        eprint(samtools_index_CMD)
        subprocess.run(samtools_index_CMD,shell=True,stdout=subprocess.PIPE)


        umi_tools_CMD="umi_tools count --extract-umi-method=tag --umi-tag=UB --cell-tag=CB --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --wide-format-cell-counts -I {} -S {}".format(sort_output,umi_tools_output)
        eprint(umi_tools_CMD)
        subprocess.run(umi_tools_CMD, shell=True,stdout=subprocess.PIPE)

        umitoolsdf=pd.read_csv(umi_tools_output,delimiter='\t',index_col=0)
        umitoolsdfT=umitoolsdf.T 
        adata=ad.AnnData(umitoolsdfT)

        adata_output=os.path.join(self.count_out_dir,'allAPA.h5ad')
        adata.write(adata_output)

        print('produce h5ad Time elapsed %.2f min' %((time.time()-ctime)/60))

        return adata









    





    













    


















