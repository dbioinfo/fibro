#!/bin/bash

#set the input variables
genome="HLictTri2"
condition="temp"
samples=$(ls rawdata/"$genome"/bams/"$condition"/ | cut -f 1 -d '.')

#setup dir structure
mkdir -p data/out_gtf/"$genome"/"$condition"
mkdir -p data/counts/"$genome"/"$condition"

#stringtie each sample
for sample_name in $samples; do

    echo $sample_name
    stringtie -o data/out_gtf/"$genome"/"$condition"/"$sample_name"_stringtie.gtf -A data/counts/"$genome"/"$condition"/"$sample_name"_gene_abundances.tsv -G rawdata/"$genome"/genome/geneAnnotation.filtered.gtf -e -v -p 80 rawdata/"$genome"/bams/"$condition"/"$sample_name".sorted.bam

done
