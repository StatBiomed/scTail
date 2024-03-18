from optparse import OptionParser,OptionGroup
from ..version import __version__
import sys
from ..utils.build_ref import get_PASref,get_gene_with_one_transcript,get_generef
from ..utils.get_counts import get_PAS_count
import os
import pandas as pd
import time 
import pyranges as pr
import torch

START_TIME = time.time()

def main():
    parser = OptionParser()
    parser.add_option('--gtf','-g',dest='gtf_file',default=None,help='The annotation gtf file for your analysing species.')
    parser.add_option('--cellbarcode',dest='cdrFile',default=None,help='The file include cell barcode which users want to keep in the downstream analysis.')
    parser.add_option('--fasta','-f',dest='fasta',default=None,help='The reference genome file') 
    parser.add_option('--bam','-b',dest='bam_file',default=None,help='The bam file of aligned from STAR or other single cell aligned software.')
    parser.add_option('--outdir','-o',dest='out_dir',default=None,help='The directory for output [default : $bam_file]') 
    parser.add_option('--chromoSize',dest='chromoSize',default=None,help='The file which includes chromosome length')
    parser.add_option('--species',dest='species',default=None,help="This indicates the species that you want to analysis. Only human and mouse are supportted. You should input human or mouse")

   
   
    group0=OptionGroup(parser,"Optional arguments")

    group0.add_option("--minCount",type="int",dest="minCount",default=50,
    help="Minimum UMI counts for TC in all cells [default: 50]")

    group0.add_option('--nproc','-p',type="int",dest='nproc',default=4,
    help='Number of subprocesses [default: 4]')

    group0.add_option('--device','-d',type="int",dest='device',default=0,
    help='If you are running this program in GPU, this is the card number assigned to use [default: 0]')

    group0.add_option('--maxReadCount',type="int",dest='maxReadCount',default=10000,
    help='For each gene, the maxmium read count kept for clustering [default: 10000]')
    
    group0.add_option('--densityFC',type="float",dest='densityFC',default=0,
    help="Minimum value for maximum density / minimum density [default: 0]")

    group0.add_option('--InnerDistance',type="float",dest='InnerDistance',default=100,
    help="The resolution of each cluster [default: 100]")





    parser.add_option_group(group0)
    (options, args) = parser.parse_args()

    
    #this means that if users do not input any argument, then direct produce help. then end.
    if len(sys.argv[1:]) == 0:
        print('Welcome to scTail v%s!\n'%(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)


    
    # bam file
    if options.bam_file is None:
        print("Error: Need --bam for aligned file.")
        sys.exit(1)

    # fasta file
    if options.fasta is None:
        print("Error: Need --fasta for reference genome")
        sys.exit(1)


    # chromosome size
    if options.chromoSize is None:
        print("Error: Need --chromoSize for chromosome size")
        sys.exit(1)


    # cellbarcode 
    if options.cdrFile is None:
        print("Error: Need --cellbarcode for cellbarcode size")
        sys.exit(1)


    if options.species is None:
        print("Error: Need --species to select the species you analysis")



    #output file 
    if options.out_dir is None:
        print("Warning: no outDir provided, we use $bamfilePath/scTail\n")
        out_dir = os.path.join(os.path.dirname(os.path.abspath(options.bam_file)),"scTail")
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)


    #gtf file  
    if options.gtf_file is None:
        print("Error: Need --gtf for annotation file.")
        sys.exit(1)
    else:
        gr = pr.read_gtf(options.gtf_file)
        grdf = gr.df
        ref_out_dir=os.path.join(str(out_dir),'ref_file')
        if not os.path.exists(ref_out_dir):
            os.mkdir(ref_out_dir)
        PASrefpath=get_PASref(grdf,ref_out_dir)
        PASdf=pd.read_csv(PASrefpath,delimiter='\t')
        generefpath=get_generef(grdf,PASdf,ref_out_dir)
        one_gene_transcriptPath=get_gene_with_one_transcript(grdf,ref_out_dir)



    bam_file=options.bam_file
    minCount=options.minCount
    fasta_file=options.fasta
    #cellBarcodePath=options.cdrFile
    n_proc=options.nproc
    maxReadCount=options.maxReadCount
    densityFC=options.densityFC
    InnerDistance=options.InnerDistance
    device=options.device
    chromoSizePath=options.chromoSize
    cellbarcodePath=options.cdrFile
    species=options.species


    # Check if GPU
    device = options.device if torch.cuda.is_available() else torch.device('cpu')
    print("[*] Selected device: ", device)


    

    getTSScount=get_PAS_count(PASrefpath,generefpath,fasta_file,bam_file,out_dir,n_proc,minCount,maxReadCount,densityFC,InnerDistance,device,chromoSizePath,cellbarcodePath,species)
    #scadata=getTSScount.produce_sclevel()
    scadata=getTSScount.assign_reads2()
    



