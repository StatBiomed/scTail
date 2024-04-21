from optparse import OptionParser,OptionGroup
from ..version import __version__
import os 
import pandas as pd
import subprocess
import sys


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)




def _get_sorted_bed(sampleListPath,outdir):
    sampledf=pd.read_csv(sampleListPath)

    allsamplels=[]
    for i in sampledf.iloc[:,0].tolist():
        df=pd.read_csv(i,delimiter='\t',header=None)
        allsamplels.append(df)

    alldf=pd.concat(allsamplels,axis=0)
    alldf.sort_values([0,1,2],inplace=True)

    output_bed_dir=os.path.join(outdir,'merge')
    os.makedirs(output_bed_dir,exist_ok=True)

    output_bed_path=os.path.join(outdir,'merge','sorted_allsample_cluster.bed')
    alldf.to_csv(output_bed_path,sep='\t',header=None,index=None,)
    print(output_bed_path)

    return output_bed_path



def _run_bedtools(input_bed_path):

    parent_dir=os.path.dirname(input_bed_path)
    output_mergebed_path=os.path.join(parent_dir,'merged_cluster.bed')


    bedtools_CMD="bedtools merge -i {} -s -d 40 -c 4 -o collapse > {}".format(input_bed_path, output_mergebed_path)
    eprint(bedtools_CMD)
    subprocess.run(bedtools_CMD, shell=True,stdout=subprocess.PIPE)


    mergedf=pd.read_csv(output_mergebed_path,delimiter='\t',header=None)
    mergedf['score']=0
    beddf=mergedf[[0,1,2,4,'score',3]]

    beddf.to_csv(output_mergebed_path,sep='\t',header=None,index=None)


    return output_mergebed_path



def main():
    parser=OptionParser()
    parser.add_option('--sampleList',dest='sampleList',default=None,help='The pathway of tsv file include the path of all samples')
    parser.add_option('--outdir','-o',dest='out_dir',default=None,help='The directory for output merge bed file [default : $bam_file]') 


    group0=OptionGroup(parser,"Optional arguments")

    group0.add_option("--maxDistance",type="int",dest="maxDistance",default=40,
    help="Maximum distance between clusters allowed for clusters to be merged. [default : 40]")


    parser.add_option_group(group0)
    (options, args) = parser.parse_args()


    #this means that if users do not input any argument, then direct produce help. then end.
    if len(sys.argv[1:]) == 0:
        print('Welcome to scTail-peakMerge v%s!\n'%(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)

    
        # bam file
    if options.sampleList is None:
        print("Error: Need --sampleList for sample list file.")
        sys.exit(1)


    if options.out_dir is None:
        print("Error: Need --outdir for output directionary")
        sys.exit(1)


    sampleList=options.sampleList
    outdir=options.out_dir


    output_bed_path=_get_sorted_bed(sampleList,outdir)
    output_mergebed_path=_run_bedtools(output_bed_path)


