library(tidyverse)
genome <- "HLratNor7"
ingtf <- read_tsv(paste0("rawdata/",genome,"/genome/geneAnnotation.gtf"), col_names=F)
trans <- read_delim(paste0("rawdata/",genome,"/genome/scaffold_localID2acc"), delim=' ')

colnames(trans) <- c('scaff','access')
tlist <- as.list(trans$access)
names(tlist) <- trans$scaff

ingtf <- ingtf %>% mutate(X1 = recode(ingtf$X1, !!!tlist, .default=X1)) %>% select(X1)

write_delim(ingtf, paste0("rawdata/",genome,"/genome/geneAnnotation.col1"), col_names=F, delim='\t', quote='none')
system("bash src/translate_ids.helper.sh", intern=F)
