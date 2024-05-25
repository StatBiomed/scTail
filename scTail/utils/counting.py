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



class DoCounting():
    def __init__(self,bamfilePath,outdir,nproc,maxReadCount,cellbarcodePath,pas_cluster_path):
    
        self.bamfilePath=bamfilePath
        self.outdir=outdir
        self.nproc=nproc
        self.maxReadCount=maxReadCount

        self.cellbarcode=pd.read_csv(cellbarcodePath,delimiter='\t')['cellbarcode'].values
        self.pas_cluster_path=pas_cluster_path 

        self.count_out_dir=os.path.join(outdir,'count')
        if not os.path.exists(self.count_out_dir):
            os.mkdir(self.count_out_dir)


        generefPath=os.path.join(outdir,'ref_file','ref_gene.tsv')
        self.generefdf=pd.read_csv(generefPath,delimiter='\t')









    def _fit_lognormal(self):
        start_time=time.time()

        #rewrite
        one_transcript_path=os.path.join(self.outdir,'ref_file','one_transcript_gene.tsv')
        one_transcriptdf=pd.read_csv(one_transcript_path,delimiter='\t')
        select_transcriptdf=one_transcriptdf.sample(n=10000,random_state=666)


        fragmentsizels=[]
        for index, row in select_transcriptdf.iterrows():

            samFile, _chrom = check_pysam_chrom(self.bamfilePath, row['Chromosome'])
            reads = fetch_reads(samFile, _chrom,  row['Start'], row['End'],  trimLen_max=100)            
            reads1_paired=reads['reads1']
            reads2_paired=reads['reads2']


            reads1_paired=[r for r in reads1_paired if r.get_tag('CB') in self.cellbarcode]
            reads2_paired=[r for r in reads2_paired if r.get_tag('CB') in self.cellbarcode]


            reads1_paired=[r for r in reads1_paired if r.get_tag('GX')==row['gene_id']]
            reads2_paired=[r for r in reads2_paired if r.get_tag('GX')==row['gene_id']]





            if row['Strand']=='+':
                reads1_paired=[r for r in reads1_paired if r.is_reverse==True]
                reads2_paired=[r for r in reads2_paired if r.is_reverse==False]
                

                reads1_pas=[r.reference_end for r in reads1_paired]
                reads2_pas=[r.reference_end for r in reads2_paired]
                fragment_size=[reads1_pas[i]-reads2_pas[i] for i in range(0,len(reads1_pas))]
                # print(fragment_size)
        
            elif row['Strand']=='-':
                reads1_paired=[r for r in reads1_paired if r.is_reverse==False]
                reads2_paired=[r for r in reads2_paired if r.is_reverse==True]


                # print(reads1_paired)
                # print(reads2_paired)

                reads1_pas=[r.reference_start for r in reads1_paired]
                reads2_pas=[r.reference_start for r in reads2_paired]
                fragment_size=[reads2_pas[i]-reads1_pas[i] for i in range(0,len(reads1_pas))]
                # print(fragment_size)

            fragmentsizels.append(fragment_size)


        # print(fragmentsizels)
        fragment_flattenls=[j for i in fragmentsizels for j in i]
        filtered_fragment=[i for i in fragment_flattenls if (i<500)&(i>0)]
        [shape_fit,location_fit,scale_fit]=stats.lognorm.fit(filtered_fragment)


        x=np.linspace(np.min(filtered_fragment),np.max(filtered_fragment))

        fig = plt.figure(figsize=(4.5, 3.5), dpi=300)
        plt.hist(filtered_fragment,bins=100,density=True)
        plt.plot(x,stats.lognorm.pdf(x,shape_fit,loc=location_fit,scale=scale_fit))
        plt.xlabel('fragment size')
        plt.ylabel('probability density')
        plt.title('fit by using log-normal distribution')

        fit_fig_output_path=os.path.join(self.count_out_dir,'lognorm_fit.pdf')
        fig.savefig(fit_fig_output_path, dpi=300, bbox_inches='tight')



        return shape_fit, location_fit, scale_fit



    def _assign_reads_one_gene(self,geneid,cluster_mapped_genedf,cutoff_likelihood,shape_fit,location_fit,scale_fit):

        selectdf=cluster_mapped_genedf[cluster_mapped_genedf['gene_id']==geneid]
        gene_chr=list(selectdf['gene_chr'])[0]
        gene_start=list(selectdf['gene_start'])[0]
        gene_end=list(selectdf['gene_end'])[0]
        strand=list(selectdf['gene_strand'])[0]
            
            
        samFile, _chrom = check_pysam_chrom(self.bamfilePath, gene_chr)
        reads = fetch_reads(samFile, _chrom, gene_start, gene_end,  trimLen_max=100)

        reads2=reads['reads2u']+reads['reads2']
        reads2=[r for r in reads2 if r.get_tag('GX')==geneid]
        reads2=[r for r in reads2 if r.get_tag('CB') in self.cellbarcode]


        if len(reads2)>self.maxReadCount:
            reads2=random.sample(reads2,10000)


        pas_output_ls=[]
        if strand=='+':
            reads2df=pd.DataFrame({'CB':[r.get_tag('CB') for r in reads2],'UMI':[r.get_tag('UB') for r in reads2],'reads2_Pos':[r.reference_end for r in reads2]})
            sorteddf=reads2df.sort_values(['CB','UMI','reads2_Pos'])
            sorteddf.groupby(['CB','UMI']).tail(n=1)

            reads2_pas=sorteddf['reads2_Pos'].tolist()


            for pas_index, pas in enumerate(reads2_pas):
                clusterid_ls=[]
                cluster_likelihood=[]
                for cluster_index,cluster in selectdf.iterrows():
                    pas_frag_size=cluster['cluster_end']-pas

                    if pas_frag_size>0:        
                        likelihood=stats.lognorm.pdf(pas_frag_size,shape_fit,loc=location_fit,scale=scale_fit)
                        

                        if likelihood>cutoff_likelihood:
                            clusterid_ls.append(cluster['cluster_id'])
                            cluster_likelihood.append(likelihood)

                try:        
                    idx_max=np.argmax(cluster_likelihood)
                    
                    belonged_to_cluster=clusterid_ls[idx_max]
                    pas_output_ls.append([pas_index,pas,belonged_to_cluster])
                except:
                    pass



        elif strand=='-':

            reads2df=pd.DataFrame({'CB':[r.get_tag('CB') for r in reads2],'UMI':[r.get_tag('UB') for r in reads2],'reads2_Pos':[r.reference_start for r in reads2]})
            sorteddf=reads2df.sort_values(['CB','UMI','reads2_Pos'])
            sorteddf.groupby(['CB','UMI']).head(n=1)
            reads2_pas=sorteddf['reads2_Pos'].tolist()


            for pas_index, pas in enumerate(reads2_pas):
                clusterid_ls=[]
                cluster_likelihood=[]
                for cluster_index,cluster in selectdf.iterrows():
                    pas_frag_size=pas-cluster['cluster_start']
                    
                    if pas_frag_size>0:        
                        likelihood=stats.lognorm.pdf(pas_frag_size,shape_fit,loc=location_fit,scale=scale_fit)
                        
                        if likelihood>cutoff_likelihood:
                            clusterid_ls.append(cluster['cluster_id'])
                            cluster_likelihood.append(likelihood)

                try:        
                    idx_max=np.argmax(cluster_likelihood)
                    
                    belonged_to_cluster=clusterid_ls[idx_max]
                    pas_output_ls.append([pas_index,pas,belonged_to_cluster])
                except:
                    pass


        reads2_assign_df=pd.DataFrame(pas_output_ls,columns=['reads2_index','reads_end','cluster_id'])

        onegene_cluster_info_ls=[]
        cluster_idls=list(selectdf['cluster_id'])
        for cluster_id in cluster_idls:
            one_cluster_reads2_assign_df=reads2_assign_df[reads2_assign_df['cluster_id']==cluster_id]
            reads2_index=one_cluster_reads2_assign_df['reads2_index'].tolist()

            cellbarcode=[reads2[index].get_tag('CB') for index in reads2_index]
            

            cellID,counts=np.unique(cellbarcode,return_counts=True)
            clusterdf=pd.DataFrame({'cell_id':cellID,cluster_id:counts})
            clusterdf.set_index('cell_id',inplace=True)
            onegene_cluster_info_ls.append(clusterdf)

        return onegene_cluster_info_ls

            


    
    def assign_reads2(self):
        
        # print(shape_fit)
        # print(location_fit)
        # print(scale_fit)
        positive_bed_path=self.pas_cluster_path
        shape_fit, location_fit, scale_fit=self._fit_lognormal()

        

        cutoff_x=np.linspace(0,500,500)
        cutoff_likelihood=np.max(stats.lognorm.pdf(cutoff_x,shape_fit,loc=location_fit,scale=scale_fit))/2



        refbeddf=self.generefdf
        refbeddf['score']=0
        refbeddf=refbeddf[['Chromosome','Start','End','gene_id','score','Strand']]

        genebedPath=os.path.join(self.outdir,'ref_file','gene_coordinate.bed')
        refbeddf.to_csv(genebedPath,sep='\t',header=None,index=None)
        
        cluster_mapped_gene_path=os.path.join(self.count_out_dir,'cluster_mapped_gene.bed')


        bedtools_CMD="bedtools intersect -a {} -b {} -wa -wb > {}".format(positive_bed_path, genebedPath, cluster_mapped_gene_path)
        eprint(bedtools_CMD)
        subprocess.run(bedtools_CMD, shell=True,stdout=subprocess.PIPE)

        cluster_mapped_genedf=pd.read_csv(cluster_mapped_gene_path,delimiter='\t',header=None)
        cluster_mapped_genedf.columns=['cluster_chr','cluster_start','cluster_end','cluster_id','cluster_score','cluster_strand','gene_chr','gene_start','gene_end','gene_id','gene_score','gene_strand']
        cluster_mapped_genedf['cluster_id']=cluster_mapped_genedf['cluster_id']+'*'+cluster_mapped_genedf['gene_id']
        cluster_mapped_genedf=cluster_mapped_genedf[(cluster_mapped_genedf['cluster_strand']==cluster_mapped_genedf['gene_strand'])&(cluster_mapped_genedf['cluster_chr']==cluster_mapped_genedf['gene_chr'])]


        unique_gene_id=set(cluster_mapped_genedf['gene_id'])
        #print(unique_gene_id)

        pool = mp.Pool(processes=self.nproc)

        all_gene_cluster_info_ls=[]
        for geneid in unique_gene_id:
            all_gene_cluster_info_ls.append(pool.apply_async(self._assign_reads_one_gene,(geneid,cluster_mapped_genedf,cutoff_likelihood,shape_fit,location_fit,scale_fit)))
        pool.close()
        pool.join()

        results=[res.get() for res in all_gene_cluster_info_ls]
        allgene_cluster_ls=[j for i in results for j in i]


        allcellbarcodels=[]
        for onedf in allgene_cluster_ls:
            allcellbarcodels.append(onedf.index.tolist())

        uniquecell=set([j for i in allcellbarcodels for j in i])

        finaldf=pd.DataFrame(index=list(uniquecell))

        for onedf in allgene_cluster_ls:
            cluster_column=onedf.columns[0]
            finaldf[cluster_column]=finaldf.index.map(onedf[cluster_column])

        finaldf.fillna(0,inplace=True)
        finaldf=finaldf.loc[:, (finaldf != 0).any(axis=0)]

        adata=ad.AnnData(finaldf)
        #print(adata.shape)
        #print(adata)
        #adata.write('/mnt/ruiyanhou/nfs_share2/three_primer/mouse_forelimb/test/mouse_forelimb.h5ad')
        vardf=pd.DataFrame(adata.var.copy())
        vardf.reset_index(inplace=True)
        vardf.columns=['cluster_id']

        vardf=vardf.merge(cluster_mapped_genedf,on='cluster_id')
        #print(cluster_mapped_genedf)
        print(vardf)



        selectgenedf=self.generefdf[['gene_id','gene_name']]
        vardf=vardf.merge(selectgenedf,on='gene_id')
        
        vardf.drop_duplicates(inplace=True)
        print(vardf)
        vardf.set_index('cluster_id',inplace=True,drop=True)
        
        #print(vardf)
        adata.var=vardf.copy()


        nonzero_indices = np.nonzero(adata.X)
        nonzero_values = adata.X[nonzero_indices]
        sparse_matrix = csr_matrix((nonzero_values, nonzero_indices), shape=adata.X.shape)
        adata.X=sparse_matrix


        adata.var['new_cluster_name']=adata.var['cluster_chr'].astype('str')+'_'+adata.var['cluster_start'].astype('str')+'_'+adata.var['cluster_end'].astype('str')
        adata.var['original_cluster_id']=adata.var.index
        adata.var.index=adata.var['new_cluster_name']+'*'+adata.var['gene_id']




        allcluster_adata_path=os.path.join(self.count_out_dir,'all_cluster.h5ad')
        adata.write(allcluster_adata_path)



        newvardf=adata.var.copy()


        twoclusterdf=newvardf[newvardf.duplicated('gene_id',keep=False)]
        twoclusteradata=adata[:,twoclusterdf.index]
        twocluster_adata_path=os.path.join(self.count_out_dir,'two_cluster.h5ad')
        twoclusteradata.write(twocluster_adata_path)



        return adata




