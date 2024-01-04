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
        self.count_out_dir=str(outdir)+'/count/'
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
        

    
    def _get_reads(self):
        # filter reads (standard: cigar and consecutive T or A)
        start_time=time.time()
        infile = pysam.AlignmentFile(self.bamfilePath, "rb")
        outfile = pysam.AlignmentFile(self.filteredbamfilePath, "wb", template=infile)

        for read in infile:
            if (read.is_reverse==False)&(read.cigartuples[0][0]==4)&(read.cigartuples[0][1]>=60)&(read.cigartuples[1][0]==0)&(read.cigartuples[1][1]>=5):
                if re.search('T{6,}',read.query_sequence):
                    outfile.write(read)

            elif (read.is_reverse==True)&(read.cigartuples[-1][0]==4)&(read.cigartuples[-1][1]>=60)&(read.cigartuples[-2][0]==0)&(read.cigartuples[-2][1]>=5):
                if re.search('A{6,}',read.query_sequence):
                    outfile.write(read)
        infile.close()
        outfile.close()

        pysam.index(self.filteredbamfilePath)
        print('filtering spend %.2f min' %((time.time()-start_time)/60))


        #remove PCR duplication 
        start_time=time.time()
        UMI_tools_CMD="umi_tools dedup --stdin {} --stdout {} --extract-umi-method=tag --umi-tag=UR --cell-tag=CR".format(self.filteredbamfilePath,self.pcr_removedPath)
        eprint(UMI_tools_CMD)
        subprocess.run(UMI_tools_CMD, shell=True,stdout=subprocess.PIPE)
        print('remove PCR duplication: %.2f min' %((time.time()-start_time)/60))

        #remove internel priming ??

        pysam.index(self.pcr_removedPath)

        #get PolyA site 
        getreadsFile=pysam.AlignmentFile(self.pcr_removedPath,'rb')
        reads_info=[]
        for read in getreadsFile.fetch(until_eof = True):
            if read.is_reverse:
                polyA_site_0_based=read.reference_end
                reads_info.append((read.reference_name,'-',polyA_site_0_based,read.get_tag('CB'),read.get_tag('GX'),read.cigarstring))
                
            else:
                polyA_site_0_based=read.reference_start 
                reads_info.append((read.reference_name,'+',polyA_site_0_based,read.get_tag('CB'),read.get_tag('GX'),read.cigarstring))

        readsdf=pd.DataFrame(reads_info,columns=['chr','strand','PAS','cellbarcode','gene_id','cigar'])
        readsdf=readsdf[readsdf['cellbarcode']!='-']
        pasdf=readsdf.groupby('gene_id')['PAS'].agg(list).reset_index()
        cellbarcodedf=readsdf.groupby('gene_id')['cellbarcode'].agg(list).reset_index()
        mergedf=pasdf.merge(cellbarcodedf,on='gene_id')
        mergedf['count']=mergedf['PAS'].apply(lambda x: len(x))
        mergedf.set_index('gene_id',inplace=True)
        finalreadsinfodf=mergedf[(mergedf['count']>self.minCount)&(mergedf['count']<self.maxReadCount)]

        reads_info_outputPath=os.path.join(self.count_out_dir,'reads_info.tsv')
        finalreadsinfodf.to_csv(reads_info_outputPath,sep='\t')

        return finalreadsinfodf 


    def _do_clustering(self,finalreadsinfodf,geneid):

        # do hierarchical cluster
        clusterModel = AgglomerativeClustering(n_clusters=None,linkage='average',distance_threshold=self.InnerDistance)

        posiarray=np.array(finalreadsinfodf['PAS'][geneid]).reshape(-1,1)
        CBarray=np.array(finalreadsinfodf['cellbarcode'][geneid]).reshape(-1,1)

        clusterModel=clusterModel.fit(posiarray)
        #print('finish clustering fit')

        labels=clusterModel.labels_
        label,count=np.unique(labels,return_counts=True)
        #print('finish label unique')

        selectlabel=label[count>=self.minCount]
        selectcount=count[count>=self.minCount]
        finallabel=list(selectlabel[np.argsort(selectcount)[::-1]])


        altTSSls=[]
        for i in range(0,len(finallabel)):
            altTSSls.append([posiarray[labels==finallabel[i]],CBarray[labels==finallabel[i]]])
        #print(altTSSls)
                       
        return altTSSls


    
    def _do_hierarchial_cluster(self):
        start_time=time.time()

        
        finalreadsinfodf=self._get_reads() 
        #print(len(readinfodict))

        altPASdict={}
        altPASfls=[]


        pool = mp.Pool(processes=self.nproc)
        for geneid in finalreadsinfodf.index:
            altPASfls.append(pool.apply_async(self._do_clustering,(finalreadsinfodf,geneid)))
        pool.close()
        pool.join()
        results=[res.get() for res in altPASfls]


        for geneidSec, reslsSec in zip(finalreadsinfodf.index,results):
            altPASdict[geneidSec]=reslsSec
        altPASdict={k: v for k, v in altPASdict.items() if v}

        polyA_output=os.path.join(self.count_out_dir,'cluster_peak.pkl')
        with open(polyA_output,'wb') as f:
            pickle.dump(altPASdict,f)

        print('do clustering Time elapsed %.2f min' %((time.time()-start_time)/60))

        return altPASdict




    def _filter_false_positive(self):
        start_time=time.time()
        altPASdict=self._do_hierarchial_cluster()

        #used as test
        with open('/mnt/ruiyanhou/nfs_share2/three_primer/eMSC/get_peak/GSM4506635/count/cluster_peak.pkl','rb') as f:
            altPASdict=pickle.load(f)

        # flatten the cluster (for some genes which include multiple cluster)
        clusterdict={}
        for i in altPASdict.keys():
            for j in range(0,len(altPASdict[i])):
                startpos=np.min(altPASdict[i][j][0])
                stoppos=np.max(altPASdict[i][j][0])
                clustername=str(i)+'*'+str(startpos)+'_'+str(stoppos)
                clusterdict[clustername]=altPASdict[i][j]

        modedict={}
        for mykeys in clusterdict.keys():
            modedict[mykeys]=stats.mode(clusterdict[mykeys][0])[0][0]

        modedf=pd.DataFrame(modedict,index=[0])
        high_freq_posdf=modedf.T
        high_freq_posdf.reset_index(inplace=True)
        print(high_freq_posdf)
        high_freq_posdf['gene_id']=high_freq_posdf['index'].str.split('*',expand=True)[0]
        high_freq_posdf.columns=['cluster_id','PAS','gene_id']

        genedf=self.generefdf[self.generefdf['gene_id'].isin(high_freq_posdf['gene_id'])]
        genedf=genedf[['Chromosome','Strand','gene_id','gene_name']]
        mergedf=high_freq_posdf.merge(genedf,on='gene_id')
        allPASpeak_mode_output=os.path.join(self.count_out_dir,'allPASpeak_mode.tsv')
        mergedf.to_csv(allPASpeak_mode_output,sep='\t',index=None)

        #run pre-trained deep learning model 
        test_set=scDataset(self.fastafilePath,allPASpeak_mode_output)
        test_loader = DataLoader(test_set, batch_size=1, shuffle=False)

        net=Net().to(self.device)
        pretrained_model_path=pathstr=str(Path(os.path.dirname(os.path.abspath(__file__))).parents[0])+'/model/model_last.pth'
        predict_result_output=os.path.join(self.count_out_dir,'predict_result.tsv')
        positivedf=test(test_loader,self.device,net,pretrained_model_path,predict_result_output)

        keep_dict = {key: value for key, value in clusterdict.items() if key in list(positivedf['PAS'])}

        pasdf=pd.DataFrame(list(keep_dict.keys()),columns=['PAS_name'])
        pasdf['gene_id']=pasdf['PAS_name'].str.split('*',expand=True)[0]


        gene_keep_dict={}
        for geneid in list(pasdf['gene_id'].unique()):
            subsetdf=pasdf[pasdf['gene_id']==geneid]
            pasls=[]
            for pasid in subsetdf['PAS_name']:
                pasls.append(keep_dict[pasid])
            
            gene_keep_dict[geneid]=pasls

        
        #should change to another folder and save the keep_dict to the folder 
        keepdict_output=os.path.join(self.count_out_dir,'keep_cluster.pkl')
        with open(keepdict_output,'wb') as f:
            pickle.dump(gene_keep_dict,f)

        print('Filtering false positive elapsed %.2f min' %((time.time()-start_time)/60))
        return gene_keep_dict



    def _do_annotation(self,inputpar):
        geneid=inputpar[0]
        altPASitemdict=inputpar[1]
        temprefdf=self.PASrefdf[self.PASrefdf['gene_id']==geneid]


        #use Hungarian algorithm to assign cluster to corresponding transcript
        cost_mtx=np.zeros((len(altPASitemdict),temprefdf.shape[0]))
        for i in range(len(altPASitemdict)):
            for j in range(temprefdf.shape[0]):
                cluster_val=altPASitemdict[i][0]

                #this cost matrix should be corrected
                position,count=np.unique(cluster_val,return_counts=True)
                mode_position=position[np.argmax(count)]
                cost_mtx[i,j]=np.absolute(np.sum(mode_position-temprefdf.iloc[j,5]))
        row_ind, col_ind = linear_sum_assignment(cost_mtx)
        transcriptls=list(temprefdf.iloc[col_ind,:]['transcript_id'])


        pasls=list(temprefdf.iloc[col_ind,:]['PAS'])


        transcriptdict={}
        for i in range(0,len(pasls)):
            if (pasls[i]>=np.min(altPASitemdict[i][0])) & (pasls[i]<=np.max(altPASitemdict[i][0])):
                name1=str(geneid)+'_'+str(transcriptls[i])
                transcriptdict[name1]=(altPASitemdict[row_ind[i]][0],altPASitemdict[row_ind[i]][1])
            else:
                newname1=str(geneid)+'_newPAS'
                transcriptdict[newname1]=(altPASitemdict[row_ind[i]][0],altPASitemdict[row_ind[i]][1])

        return transcriptdict



    
    def _TSS_annotation(self):
        start_time=time.time()

        keepdict=self._filter_false_positive()
        keepIDls=list(keepdict.keys())
        
        inputpar=[]
        for i in keepIDls:
            inputpar.append((i,keepdict[i]))

        pool = mp.Pool(processes=self.nproc)
        with mp.Pool(self.nproc) as pool:
            transcriptdictls=pool.map_async(self._do_annotation,inputpar).get()


        extendls=[]
        for d in transcriptdictls:
            extendls.extend(list(d.items()))


        
        ### organize the output result
        d={'transcript_id':[transcript[0] for transcript in extendls],'PAS_start':[np.min(transcript[1][0]) for transcript in extendls],
        'PAS_end':[np.max(transcript[1][0]) for transcript in extendls]}

        regiondf=pd.DataFrame(d)
        print('Doing annotation elapsed %.2f min' %((time.time()-start_time)/60))

        return extendls,regiondf


    
    def produce_sclevel(self):
        ctime=time.time()
        extendls,regiondf=self._TSS_annotation()

        cellIDls=[]
        for i in range(0,len(extendls)):
            cellID=np.unique(extendls[i][1][1])
            cellIDls.append(list(cellID))
        cellIDset = set([item for sublist in cellIDls for item in sublist])
        finaldf=pd.DataFrame(index=list(cellIDset))



        for i in range(0,len(extendls)):
            transcriptid=extendls[i][0]       
            cellID,count=np.unique(extendls[i][1][1],return_counts=True)
            transcriptdf=pd.DataFrame({'cell_id':cellID,transcriptid:count})
            transcriptdf.set_index('cell_id',inplace=True)
            finaldf[transcriptid]=finaldf.index.map(transcriptdf[transcriptid])


        finaldf.fillna(0,inplace=True)
        adata=ad.AnnData(finaldf)
        vardf=pd.DataFrame(adata.var.copy())
        vardf.reset_index(inplace=True)
        vardf.columns=['transcript_id']
        vardf=vardf.join(regiondf.set_index('transcript_id'), on='transcript_id')
        vardf['gene_id']=vardf['transcript_id'].str.split('_',expand=True)[0]
        vardf=vardf.merge(self.generefdf,on='gene_id')
        vardf.set_index('transcript_id',drop=True,inplace=True)

        adata.var=vardf
        sc_output_h5ad=os.path.join(self.count_out_dir,'scPAS_count_all.h5ad')
        adata.write(sc_output_h5ad)

        #filter according to user' defined distance
        newdf=adata.var.copy()
        newdf.reset_index(inplace=True)
        selectedf=newdf[newdf.duplicated('gene_id',keep=False)]  #get data frame which includes two transcript for one gene
        geneID=selectedf['gene_id'].unique()

        keepdfls=[]
        for i in geneID:
            tempdf=selectedf[selectedf['gene_id']==i]

            tempdf=tempdf.sort_values('transcript_id',ascending=False)
            tempdf['diff']=tempdf['PAS_start'].diff()
            keepdf=tempdf[tempdf['diff'].isna()|tempdf['diff'].abs().ge(self.clusterDistance)]    #want to get TSS whose cluster distance is more than user defined.
            #keepdf=keepdf.iloc[:2,:]
            keepdfls.append(keepdf) 

        #print(keepdfls)


        allkeepdf=reduce(lambda x,y:pd.concat([x,y]),keepdfls)
        finaltwodf=allkeepdf[allkeepdf.duplicated('gene_id',keep=False)] 
        finaltwoadata=adata[:,adata.var.index.isin(finaltwodf['transcript_id'])]  

        sc_output_h5ad=os.path.join(self.count_out_dir+'scPAS_count_two.h5ad')
        finaltwoadata.write(sc_output_h5ad)

        print('produce h5ad Time elapsed %.2f min' %((time.time()-ctime)/60))

        return adata






    





    













    


















