#!/bin/bash

#set the input variables
genome="HLictTri2"
condition="glucose"
samples=$(ls rawdata/"$genome"/fastq/ | cut -f 1 -d '_' | sort | uniq)

#build an index for the assembly
#hisat2-build -p 80 rawdata/"$genome"/genome/"$genome".fa rawdata/"$genome"/genome/"$genome" 

mkdir -p rawdata/"$genome"/bams/"$condition"

#align each sample
for sample_name in $samples; do

    echo $sample_name
    hisat2 -p 80 --dta -x rawdata/"$genome"/genome/"$genome" -1 rawdata/"$genome"/fastq/"$sample_name"_R1_001.fastq.gz -2 rawdata/"$genome"/fastq/"$sample_name"_R2_001.fastq.gz | samtools sort -@ 80 -o rawdata/"$genome"/bams/"$condition"/"$sample_name".sorted.bam


done
