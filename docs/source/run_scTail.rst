==============
Run scTail
==============

scTail includes three stpes : identify PAS for individual samples, merge PASs across multiple samples, quantify PAS for individual cells  

The input files include:

* filtered alignment file (bam file)
* annotation file (gtf file)
* reference genome file (fasta file)
* cell barcode list file (tsv file)
* chrom size (tsv file)

The output files include:

* paraclu_input.tsv : a file includes 4 columns: chromosome, strand, position, reads_count; this file is a middle file and will input to paraclu.
* paraclu_output.tsv : a file contains 8 columns: the sequence name, the strand, the first position in the cluster, the last position in the cluster, the number of positions with data in the cluster, the sum of the data values in the cluster, the cluster's "minimum density" and the cluster's "maximum density".    
* input_to_DP.tsv : a file include 4 columns: chromosome, strand, PAS and cluster_id. This file was made to input to our deep learning neural network.
* predict_result.tsv : a file contains three columns: cluster_id; predicted label by CNN model; probability to be positive sample. 
* positive_result.bed : A dataframe contains six columns: chromosome, cluster_start, cluster_end, cluster_name, score and strand. Samples in this file are positive samples. 
* lognorm_fit.pdf : a file shows the distribution of fragment size and the result by using log normal distribution to fit
* cluster_mapped_gene.bed : a file contains 11 columns: chromosome of PAS cluster, start of PAS cluster, end of PAS cluster, PAS cluster id, score of PAS cluster, strand of PAS cluster, chromosome of gene, start of gene, end of gene, gene id, score
* all_cluster.h5ad: final result, including matrix of cell by all PAS. In most situations, you just need focus this file.
* two_cluster.h5ad: including matrix of cell by alternative PAS (with two or more than PAS).


Here is a quick test file. You can check it.
  
Download test file
===================

You can download test file from figshare_.

.. _figshare: https://doi.org/10.6084/m9.figshare.25902508.v1

Here, you can download some large file include genome.fa, annotation gtf, a bam file and so on. Please note that you should also download the index file (.bam.bai) for the bam file or index by youself.

Alternatively, you can also download the reference genome fasta file from Ensembl or Genecode. But please make you use the same genome.fa and annotation gtf file as the process building STAR index.

Run scTail
=============

Here are three stpes in scTail : **scTail-callPeak** , **scTail-peakMerge** and **scTail-count**.

* scTail-callPeak : identify PAS signal for individual bam file 
* scTail-peakMerge : merge multiple PAS signal from different samples.
* scTail-count : quantify PAS signal for individual cells. 



You can run scTail-callPeak by using test file according to the following code.

.. code-block:: bash

   #!/bin/bash
   gtfFile=$download/gencode.v44.annotation.gtf
   fastaFile=$download/GRCh38.primary_assembly.genome.fa
   bamFile=$download/PBMC866_GX_CB_UB_filtered_test.bam
   cellbarcodeFile=$download/barcodes.tsv
   chromosomeSize=$download/hg38.chromsize
   outputfile=your_faverate_output

   scTail-callPeak -b $bamFile --gtf $gtfFile --cellbarcode $cellbarcodeFile -f $fastaFile --species human --chromoSize $chromosomeSize -o $outputfile --minCount 50 -p 20

   
 

Before running scTail-peakMerge, you should create a positivebed_list.tsv. This file looks like this.

.. code-block:: bash

    $your_faverate_output/count/positive_result.bed
    $download/onesample_positive.bed


Afterwards, you can run scTail-peakMerge by using test file according to the following code.

.. code-block:: bash

   #!/bin/bash
   positivebed_list=get_by_yourself
   outputfile=same_with_the_output_of_scTail-callpeak

   scTail-peakMerge --sampleList $positivebed_list.tsv -o $outputfile 


You can run scTail-count by using test file according to the following code. 

.. code-block:: bash

   #!/bin/bash 
   cellbarcodeFile=$download/barcodes.tsv
   bamFile=$download/PBMC866_GX_CB_UB_filtered_test.bam
   mergebedFile=$outputfile/merge/merged_cluster.bed
   outputfile=same_with_the_output_of_scTail-callpeak
   
   scTail-count --cellbarcode $cellbarcodeFile --bam $bamFile --outdir $outputfile --PAScluster $mergebedFile 



Options
========


There are more parameters for setting (``scTail-callPeak -h`` always give the version
you are using):


.. code-block:: html

   Usage: scTail-callPeak [options]

   Options:
        -h, --help            show this help message and exit
        -g GTF_FILE, --gtf=GTF_FILE
                        The annotation gtf file for your analysing species.
        --cellbarcode=CELLBARCODE
                        The file include cell barcode which users want to keep
                        in the downstream analysis.
        -f FASTA, --fasta=FASTA
                        The reference genome file
        -b BAM_FILE, --bam=BAM_FILE
                        The bam file of aligned from STAR or other single cell
                        aligned software.
        -o OUT_DIR, --outdir=OUT_DIR
                        The directory for output [default : $bam_file]
        --chromoSize=CHROMOSIZE
                        The file which includes chromosome length
        --species=SPECIES     This indicates the species that you want to analysis.
                        Only human and mouse are supportted. You should input
                        human or mouse

   Optional arguments:
        --minCount=MINCOUNT
                        Minimum UMI counts for one cluster in all cells
                        [default: 50]
        -p NPROC, --nproc=NPROC
                        Number of subprocesses [default: 4]
        -d DEVICE, --device=DEVICE
                        If your server has the GPU, then the default card 0
                        will be used. If your server did not have the GPU,
                        then cpu will be used.
        --maxReadCount=MAXREADCOUNT
                        For each gene, the maxmium read count kept for
                        clustering [default: 10000]
        --densityFC=DENSITYFC
                        Minimum value for maximum density / minimum density
                        [default: 0]
        --InnerDistance=INNERDISTANCE
                        The resolution of each cluster [default: 100]


There are more parameters for setting (``scTail-peakMerge -h`` always give the version
you are using):


.. code-block:: html

   Usage: scTail-peakMerge [options]

   Options:
        -h, --help            show this help message and exit
        --sampleList=SAMPLELIST
                        The pathway of tsv file include the path of all
                        samples
        -o OUT_DIR, --outdir=OUT_DIR
                        The directory for output merge bed file [default :
                        $bam_file]

   Optional arguments:
        --maxDistance=MAXDISTANCE
                        Maximum distance between clusters allowed for clusters
                        to be merged. [default : 40]



There are more parameters for setting (``scTail-count -h`` always give the version
you are using):


.. code-block:: html

   Usage: scTail-count [options]

   Options:
        -h, --help            show this help message and exit
        --cellbarcode=CELLBARCODE
                        The file include cell barcode which users want to keep
                        in the downstream analysis.
        -b BAM_FILE, --bam=BAM_FILE
                        The bam file of aligned from STAR or other single cell
                        aligned software.
        -o OUT_DIR, --outdir=OUT_DIR
                        The directory for output [default : $bam_file]
        --PAScluster=PASCLUSTER
                        The bed file of PAS cluster

   Optional arguments:
        -p NPROC, --nproc=NPROC
                        Number of subprocesses [default: 4]
        --maxReadCount=MAXREADCOUNT
                        For each gene, the maxmium read count kept for
                        clustering [default: 10000]



