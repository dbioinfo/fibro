library(rtracklayer)
library(plyranges)
library(stringr)
library(readr)

conditions <- c('glucose','hypoxia','temp')

for (ic in conditions){

    fnames <- list.files(paste0('data/out_gtf/hg38/',ic,'/'))

    for (fname in fnames){ 

        gtf <- import(paste0('data/out_gtf/hg38/',ic,'/', fname))

        out <- gtf %>% 
            filter(type=='transcript') %>% as.data.frame() %>% 
            mutate("Gene ID" = paste(transcript_id, ref_gene_name, sep='.'), Reference=1) %>% 
            rename("gene_id"="Gene Name", "start"="Start", "end"="End", "cov"="Coverage", "strand"="Strand") %>% 
            select(c("Gene ID", "Gene Name", "Reference", "Strand", "Start", "End", "Coverage", "FPKM", "TPM")) 


        outname <- paste0('data/counts/hg38/',ic,'/', str_split(fname, pattern='\\.')[[1]][1], '_gene_abundances.tsv')
        write_tsv(out, outname)

    }

}
