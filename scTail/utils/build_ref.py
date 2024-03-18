import pandas as pd
import pyranges as pr
import numpy as np
from functools import reduce
import os


def get_PASref(grdf,out_dir):
    transcriptdf=grdf[grdf['Feature']=='transcript']
    transcriptdf=transcriptdf[['transcript_id','gene_id','gene_name','Chromosome','Start','End','Strand']]
    transcriptdf['PAS']=np.where(transcriptdf['Strand']=='+',transcriptdf['End'],transcriptdf['Start'])
    transcriptdf.drop_duplicates(['Chromosome','PAS'],inplace=True)
    transcriptdf=transcriptdf[['transcript_id','gene_id','gene_name','Chromosome','Strand','PAS']]
    transcriptdf.dropna(subset=['Chromosome'],axis=0,inplace=True)

    if transcriptdf['Chromosome'].str.contains('chr').any():
        transcriptdf['Chromosome']=transcriptdf['Chromosome']
    else:
        transcriptdf['Chromosome']='chr'+transcriptdf['Chromosome']

    PASoutput_path=os.path.join(out_dir,'ref_PAS.tsv')
    transcriptdf.to_csv(PASoutput_path,index=None,sep='\t')
    return PASoutput_path


def get_generef(grdf,transcriptdf,out_dir):
    genedf=grdf[grdf['Feature']=='gene']
    genedf=genedf[['Chromosome','Feature','Start','End','Strand','gene_id','gene_name']]

    if genedf['Chromosome'].str.contains('chr').any():
        genedf['Chromosome']=genedf['Chromosome']
    else:
        genedf['Chromosome']='chr'+genedf['Chromosome']

    genedf=genedf[genedf['gene_id'].isin(pd.unique(transcriptdf['gene_id']))]
    genedf.dropna(subset=['Chromosome'],axis=0,inplace=True)

    geneoutput_path=os.path.join(out_dir,'ref_gene.tsv')
    genedf.to_csv(geneoutput_path,index=None,sep='\t')
    return geneoutput_path


def get_gene_with_one_transcript(grdf,out_dir):
    transcriptdf=grdf[grdf['Feature']=='transcript']
    transcriptdf=transcriptdf[['Chromosome','Start','End','Strand','gene_id','transcript_id','gene_name','transcript_name']]
    countdf=transcriptdf.groupby('gene_id').size().reset_index(name='count')
    onetranscriptdf=countdf[countdf['count']==1]
    onetranscriptbed=transcriptdf[transcriptdf['gene_id'].isin(onetranscriptdf['gene_id'])]

    onetranscript_path=os.path.join(out_dir,'one_transcript_gene.tsv')
    onetranscriptbed.to_csv(onetranscript_path,sep='\t',index=None)
    return onetranscript_path


