#!/bin/bash

scTail \
    -b /mnt/ruiyanhou/nfs_share2/three_primer/PBMC/preprocess/STAR_output/PBMC866_GX_CB_UB_filtered_test.bam \
    --gtf /mnt/ruiyanhou/nfs_share2/annotation/annotation_from_Gencode/gencode.v44.annotation.gtf \
    --cellbarcode /mnt/ruiyanhou/nfs_share2/three_primer/PBMC/preprocess/STAR_output/PBMC866Solo.out/Gene/filtered/barcodes.tsv \
    -f /mnt/ruiyanhou/nfs_share2/annotation/annotation_from_Gencode/GRCh38.primary_assembly.genome.fa \
    -d 1 \
    --chromoSize /mnt/ruiyanhou/nfs_share2/annotation/annotation_from_Gencode/hg38.chromsize \
    -o /mnt/ruiyanhou/nfs_share2/three_primer/PBMC/run_scTail/PBMC866_final_extract_reads_test \
    --minCount 50 \
    -p 40


