from optparse import OptionParser,OptionGroup
from ..version import __version__
import sys
from ..utils.counting import DoCounting
import os
import pandas as pd
import time 


START_TIME = time.time()

def main():
    parser = OptionParser()
    parser.add_option('--cellbarcode',dest='cellbarcode',default=None,help='The file include cell barcode which users want to keep in the downstream analysis.')
    parser.add_option('--bam','-b',dest='bam_file',default=None,help='The bam file of aligned from STAR or other single cell aligned software.')
    parser.add_option('--outdir','-o',dest='out_dir',default=None,help='The directory for output [default : $bam_file]') 
    parser.add_option('--PAScluster',dest='PAScluster',default=None,help='The bed file of PAS cluster')

   
   
    group0=OptionGroup(parser,"Optional arguments")



    group0.add_option('--nproc','-p',type="int",dest='nproc',default=4,
    help='Number of subprocesses [default: 4]')


    group0.add_option('--maxReadCount',type="int",dest='maxReadCount',default=10000,
    help='For each gene, the maxmium read count kept for clustering [default: 10000]')
    



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




    # chromosome size
    if options.PAScluster is None:
        print("Error: Need --PAScluster for PAS cluster bed file")
        sys.exit(1)


    # cellbarcode
    if options.cellbarcode is None:
        print("Error: Need --cellbarcode for cellbarcode tsv file")
        sys.exit(1)






    #output file 
    if options.out_dir is None:
        print("Warning: no outDir provided, we use $bamfilePath/scTail\n")
        out_dir = os.path.join(os.path.dirname(os.path.abspath(options.bam_file)),"scTail")
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)





    bam_file=options.bam_file
    cellBarcodePath=options.cellbarcode
    n_proc=options.nproc
    maxReadCount=options.maxReadCount
    PAScluster=options.PAScluster

    

    getPAScount=DoCounting(bam_file,out_dir,n_proc,maxReadCount,cellBarcodePath,PAScluster)
    scadata=getPAScount.assign_reads2()
    



