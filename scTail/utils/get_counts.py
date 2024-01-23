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
    def __init__(self,PASrefPath,generefPath,fastafilePath,bamfilePath,outdir,nproc,minCount,maxReadCount,clusterDistance,InnerDistance,device):
    
        self.PASrefdf=pd.read_csv(PASrefPath,delimiter='\t')
        self.generefdf=pd.read_csv(generefPath,delimiter='\t')
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





    def _do_preprocess(self):
        start_time=time.time()

        filteredbamfilePath=self.filteredbamfilePath
        infile = pysam.AlignmentFile(self.bamfilePath, "rb")
        outfile = pysam.AlignmentFile(filteredbamfilePath, "wb", template=infile)
        for read in infile:
            if read.get_tag('GX')!='-':
                outfile.write(read)
        infile.close()
        outfile.close()

        pysam.index(self.filteredbamfilePath)
        print('filtering reads elapsed %.2f min' %((time.time()-start_time)/60))

        return filteredbamfilePath




    def _get_reads(self):

        start_time=time.time()
        filteredbamfilePath=self._do_preprocess()


        start_time=time.time()
        UMI_tools_CMD="umi_tools dedup --stdin {} --stdout {} --extract-umi-method=tag --umi-tag=UR --cell-tag=CR".format(filteredbamfilePath,self.pcr_removedPath)
        eprint(UMI_tools_CMD)
        subprocess.run(UMI_tools_CMD, shell=True,stdout=subprocess.PIPE)
        print('remove PCR duplication: %.2f min' %((time.time()-start_time)/60))

        samtools_index_CMD="samtools index -@ {} {}".format(self.nproc,self.pcr_removedPath)
        eprint(samtools_index_CMD)
        subprocess.run(samtools_index_CMD,shell=True,stdout=subprocess.PIPE)



        getreadsFile=pysam.AlignmentFile(self.pcr_removedPath,'rb')
        geneidls=[]
        for read in getreadsFile.fetch(until_eof = True):
            geneid=read.get_tag('GX')
            geneidls.append(geneid)

        geneiddf=pd.DataFrame(geneidls,columns=['gene_id'])
        geneid_uniqdf=geneiddf.drop_duplicates('gene_id')

        mergedf=geneid_uniqdf.merge(self.generefdf,on='gene_id')
        mergedf.set_index('gene_id',inplace=True)
        #print(mergedf)


        geneIDls=[]
        reads1_POS_ls=[]
        cellbarcodels=[]
        for geneid in mergedf.index:
            geneIDls.append(geneid)

            #print(geneid)
            mysamFile=self.pcr_removedPath
            #print(mysamFile)
            samFile, _chrom = check_pysam_chrom(mysamFile, str(mergedf.loc[geneid]['Chromosome']))
            
            reads = fetch_reads(samFile, _chrom,  mergedf.loc[geneid]['Start'] , mergedf.loc[geneid]['End'],  trimLen_max=100)
            reads1_umi = reads["reads1u"]
            reads1_umi=[r for r in reads1_umi if r.get_tag('GX')==geneid]

            if mergedf.loc[geneid]['Strand']=='+':
                reads1_umi=[r for r in reads1_umi if r.is_reverse==False]
                reads1_POS=[r.reference_end for r in reads1_umi]
                cellbarcode=[r.get_tag('CB') for r in reads1_umi]
                reads1_POS_ls.append(reads1_POS)
                cellbarcodels.append(cellbarcode)
                
            elif mergedf.loc[geneid]['Strand']=='-':
                reads1_umi=[r for r in reads1_umi if r.is_reverse==True]
                reads1_POS=[r.reference_start for r in reads1_umi]
                cellbarcode=[r.get_tag('CB') for r in reads1_umi]
                reads1_POS_ls.append(reads1_POS)
                cellbarcodels.append(cellbarcode)

        readsinfodict={}
        readsinfodict['gene_id']=geneIDls
        readsinfodict['reads1']=reads1_POS_ls
        readsinfodict['cellbarcode']=cellbarcodels

        readsinfodf=pd.DataFrame(readsinfodict)
        readsinfodf['count']=readsinfodf['reads1'].apply(lambda x: len(x))


        #get the PAS of genes which only have one transcript
        trained_pas_path=str(Path(os.path.dirname(os.path.abspath(__file__))).parents[0])+'/model/human_hg38_gene_withone_PAS.tsv'
        trainedpasdf=pd.read_csv(trained_pas_path,delimiter='\t')

        trainedpasdf['gene_id']=trainedpasdf['gene_id'].str.split('.',expand=True)[0]
        readsinfodf['gene_id']=readsinfodf['gene_id'].str.split('.',expand=True)[0]
        print(trainedpasdf)
        trainedfulldf=readsinfodf.merge(trainedpasdf,on='gene_id')
        trainedfulldf.set_index('gene_id',inplace=True)
        print(trainedfulldf)

        #get the distance between reads and PAS for these genes 
        distancedict={}
        for geneid in trainedfulldf.index:
            distancels=[]
            for i in range(0,len(trainedfulldf.loc[geneid]['reads1'])):
                if trainedfulldf.loc[geneid]['strand']=='+':
                    distance=trainedfulldf.loc[geneid]['PAS']-trainedfulldf.loc[geneid]['reads1'][i]
                elif trainedfulldf.loc[geneid]['strand']=='-':
                    distance=trainedfulldf.loc[geneid]['reads1'][i]-trainedfulldf.loc[geneid]['PAS']
                distancels.append(distance)
            distancedict[geneid]=distancels

        #print(distancedict)

        trainedfulldf['distance']=distancedict
        train_dataset_Path=os.path.join(self.count_out_dir,'used_for_train_dataset.tsv')
        trainedfulldf.to_csv(train_dataset_Path,sep='\t',index=None)


        # use gaussian kernel function to estimate the distribution density
        flattenls=[j for i in trainedfulldf['distance'] for j in i ]
        distancenp=np.array(flattenls).reshape([-1,1])
        selectnp=distancenp[distancenp<2000].reshape(-1,1)
        kde = KernelDensity(kernel='gaussian', bandwidth=1).fit(selectnp)

        x_range = np.linspace(0,2000, num=2000)
        log_density = kde.score_samples(x_range.reshape(-1,1))

        fragmentLen=np.argmax(log_density)
        print("The predicted fragment length is %i"%fragmentLen)

        corrected_reads1dict={}
        for geneid in readsinfodf.index:
            reads1=readsinfodf.loc[geneid]['reads1']
            corrected_reads1dict[geneid]=[ele+fragmentLen for ele in reads1]
        readsinfodf['corrected_reads1']=corrected_reads1dict

        readsinfodf=readsinfodf[readsinfodf['count']>self.minCount]
        reads_info_outputPath=os.path.join(self.count_out_dir,'corrected_reads_info.tsv')
        readsinfodf.to_csv(reads_info_outputPath,sep='\t')

        print('getting reads elapsed %.2f min' %((time.time()-start_time)/60))

        return readsinfodf 




    
    def _do_cluster(self):
        start_time=time.time()

        
        readsinfodf=self._get_reads() 

        allrowls=[]
        for index, row in readsinfodf.iterrows():
            gene_id = row['gene_id']
            corrected_reads = row['corrected_reads1']
            
            # statistic the element in the list 
            counts = Counter(corrected_reads)
            
            # extract the unique value and corresponding count
            unique_values = list(counts.keys())
            counts = list(counts.values())

            onerowdf=pd.DataFrame({'gene_id': gene_id, 'unique_values': unique_values, 'count': counts})
            allrowls.append(onerowdf)
            # add all results to the new dataframe

        paracluinputdf=pd.concat(allrowls)
        paracluinputdf=paracluinputdf.explode(['unique_values','count'])
        finalparcluinputdf=paracluinputdf.merge(self.generefdf,on='gene_id')
        finalparcluinputdf=finalparcluinputdf[['Chromosome','Strand','unique_values','count']]
        finalparcluinputdf.sort_values(['Chromosome','Strand','unique_values'],inplace=True)

        paraclu_inputPath=os.path.join(self.count_out_dir,'paraclu_input.tsv')
        finalparcluinputdf.to_csv(paraclu_inputPath,sep='\t',header=None,index=None)

        paraclu_outputPath=os.path.join(self.count_out_dir,'paraclu_output.tsv')

        paraclu_CMD="paraclu {} {} | paraclu-cut > {}".format(self.minCount, paraclu_inputPath, paraclu_outputPath)
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
        test_set=scDataset(self.fastafilePath,input_to_DP)
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









    





    













    


















