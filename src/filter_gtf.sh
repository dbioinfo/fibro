#!/bin/bash

genome="HLictTri2"
ingtf="rawdata/$genome/genome/geneAnnotation.gtf"
orthos="rawdata/$genome/genome/orthologsClassification.tsv"
outgtf="rawdata/$genome/genome/geneAnnotation.filtered.gtf"

#filter the ortho file for 1-1 genes only (5th column), return IDs (4th column) to a temp file
awk -v q="\"" -F'\t' '$5 == "one2one" {print q $4 q}' "$orthos" > "temp.txt"

#search gtf for lines with matching (literal) IDs, output to new file
grep -F -f "temp.txt" "$ingtf" > "$outgtf"

#cleanup temps
rm temp.txt
