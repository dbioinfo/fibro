genome="HLratNor7"
cut -f2- rawdata/"$genome"/genome/geneAnnotation.gtf > rawdata/"$genome"/genome/geneAnnotation.col2
paste rawdata/"$genome"/genome/geneAnnotation.col1 rawdata/"$genome"/genome/geneAnnotation.col2 > rawdata/"$genome"/genome/geneAnnotation.translated.gtf
rm rawdata/"$genome"/genome/geneAnnotation.col1
rm rawdata/"$genome"/genome/geneAnnotation.col2
mv rawdata/"$genome"/genome/geneAnnotation.gtf rawdata/"$genome"/genome/geneAnnotation.raw.gtf
mv rawdata/"$genome"/genome/geneAnnotation.translated.gtf rawdata/"$genome"/genome/geneAnnotation.gtf

