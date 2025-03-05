library(tidyverse)
library(DESeq2)
library(RColorBrewer)

graph_single_gene <- function(gene, icounts, meta){
  
  #df <- icounts[gene,] %>% 
  #  pivot_longer(cols=everything(), names_to='sample', values_to = 'counts') %>% 
  #  left_join(meta, by='sample')
  
  df <- data.frame(counts=icounts[gene,], sample=names(icounts[gene,])) %>% 
    left_join(meta, by='sample') %>% group_by(individual)
  
  ggplot(df)+
    geom_jitter(mapping=aes(y=counts, x=individual, color=as.factor(treatment), shape=as.factor(treatment)), 
                width=0.1, height=0, size=3)+ 
    ylab('Normalized Counts')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0(gene, ' Normalized Counts'))+
    theme_bw()
  
}

graph_library_sizes <- function(icounts, meta){
  tmp <- data.frame(library_size= colSums(icounts[,2:length(icounts)], na.rm = T))
  tmp$sample <- rownames(tmp)
  tmp <- tmp %>% left_join(., meta, by="sample") %>% 
    group_by(individual, treatment) %>% 
    mutate(treatment = factor(treatment, levels=c("8","2.5","30","6","24","32","41")))
  ggplot(tmp)+
    geom_bar(mapping=aes(x=sample, y=library_size, fill=treatment), stat = 'identity')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_fill_manual(values =  c("#bdbdbd", brewer.pal(length(unique(tmp$treatment)) - 1 , "Paired")) )
  
}

graph_PCA <- function(icounts, imeta){
  PCA <- princomp(counts)
  PCA <- as.data.frame(PCA$loadings[,1:2])
  PCA$sample <- rownames(PCA)
  
  df <- left_join(imeta, PCA, by='sample') %>% mutate(treatment=factor(treatment))
  ggplot(df)+
    geom_text(mapping=aes(x=Comp.1, y=Comp.2, label=sample, color=treatment),
              size=5)+
    theme_bw()
  
  
}


setwd("~/Postdoc/Hindle/dylan/fibro")

#list all count matrices
fnames <- list.files('data/counts/full_transcript_counts/', full.names = T)

#read in
allcounts <- imap(fnames, read.csv)
names(allcounts) <- str_split_i(gsub("\\.csv", "", basename(fnames)), '_', 1)

#OPTIONALLY collapse transcripts to genes
allcounts <- allcounts %>% map( ~ mutate(.x, Transcript = str_split_i(Transcript, '\\.', 2)) %>% 
                                  group_by(Transcript) %>% 
                                  summarise_all(sum) %>% 
                                  ungroup())

#format or remove bad samples
#allcounts$fruitbat <- allcounts$fruitbat %>% rename("RA142p30mM"="RA142p30mMpRNAp1.18.24") #some weird sample names
#allcounts$squirrel <- allcounts$squirrel %>% rename("FIT60p24H" = "FIT60p24h", 
#                                                    "FIT60p6H" = "FIT60p6h") 
#allcounts$rat <- allcounts$rat %>% rename("RN30.2.5mM" = "RN30.2.5.mM")

allcounts$fruitbat <- allcounts$fruitbat %>% rename("RA142p30mMpRNAp1.18.24"="RA142p30mM") #stupid rename function
allcounts$squirrel <- allcounts$squirrel %>% rename("FIT60p24h"="FIT60p24H" , 
                                                    "FIT60p6h"= "FIT60p6H") 
allcounts$rat <- allcounts$rat %>% rename("RN30.2.5.mM"="RN30.2.5mM" )



allcounts$human <- allcounts$human %>% select(-c(HS967.2.5mM)) #bad samples
allcounts$lbbat <- allcounts$lbbat %>% select(-c(PG554.8mM, PG556.6H)) 
allcounts$rat <- allcounts$rat %>% select(-c(RN28.8mM, RN28.32C, RN28.41C))
allcounts$squirrel <- allcounts$squirrel %>% select(-c(OK29.6H))

#make meta from sample names
cnames <- c()
for(i in 1:length(allcounts)){
  cnames <- c(cnames, colnames(allcounts[[i]])[2:length(allcounts[[i]])])
}
meta <- tibble(sample=cnames) %>% 
  mutate(individual = str_split_i(sample, "\\.|p", 1)) %>% 
  mutate(exp_group = case_when(
    str_sub(sample, -1, -1) == 'C' ~ "temperature",
    str_sub(sample, -1, -1) == 'H' ~ "hypoxia",
    str_sub(sample, -1, -1) == 'M' ~ "glucose",
    .default = NA
  )) %>% 
  mutate(treatment = case_when(
    str_sub(sample, -3, -1) == '5mM' ~ 2.5,
    str_sub(sample, -3, -1) == '8mM' ~ 8,
    str_sub(sample, -4, -1) == '30mM' ~ 30,
    str_sub(sample, -2, -1) == '6H' ~ 6,
    str_sub(sample, -3, -1) == '24H' ~ 24,
    str_sub(sample, -3, -1) == '32C' ~ 32,
    str_sub(sample, -3, -1) == '41C' ~ 41,
    .default = NA
  ))

treats <- c('glucose', 'hypoxia', 'temperature')
treatlevs <- list('glucose'=c(8,2.5,30), 'hypoxia'=c(8,6,24), 'temperature'=c(8, 32, 41))
#do *only minimal* pre-processing, start with RAW COUNTS, don't forget 8mM glucose is control
i <- 1
j <- 'temperature'
for (i in 1:length(allcounts)){ #iter species
  for (j in treats){ #iter treatment
    if (nrow(meta %>% filter(sample %in% colnames(allcounts[[i]]), exp_group==j)) >0){ #make sure we have data for combo
      print(paste0(names(allcounts[i]), " ", j))
      counts <- as.data.frame(allcounts[[i]]) #format subset of counts for limma
      rownames(counts) <- counts$Transcript 
      counts <- counts %>% select(-Transcript)
      
      if (j=='glucose'){imeta <- meta %>% filter(sample %in% colnames(counts), exp_group==j) #subset meta
      } else { imeta <- meta %>% filter(sample %in% colnames(counts), exp_group==j) 
      imeta <- rbind(imeta, meta %>% filter(sample %in% colnames(counts), treatment==8))} #get 8mM control
      
      counts <- counts[,match(imeta$sample, colnames(counts))] #make absolute sure the cols and rows (samples) match order
      counts[is.na(counts)] <- 0 #no NA (why would NA even be possible when 0 is?)
      
      imeta$treatment <- factor(as.character(imeta %>% pull(treatment)),  levels =  treatlevs[[j]])
      mincount <- 5
      gidx <- rowSums( counts >= mincount ) >= 3
      dds <- DESeqDataSetFromMatrix(countData = counts[gidx,],
                                    colData = imeta,
                                    design = ~individual + treatment)
      
      dds <- estimateSizeFactors(dds)
      dds <- DESeq(dds)
      
      #write_csv(results(dds, contrast = c("treatment", treatlevs[[j]][2],"8"))%>% 
      #            as.data.frame() %>% arrange(padj) %>% rownames_to_column('gene'),
      #          paste0('data/DESeq/', names(allcounts[i]), "_", j, "_",treatlevs[[j]][2],"v8_DEG.csv") ) 
      #write_csv(results(dds, contrast = c("treatment", treatlevs[[j]][3],"8"))%>% 
      #            as.data.frame() %>% arrange(padj) %>% rownames_to_column('gene'),
      #          paste0('data/DESeq/', names(allcounts[i]), "_", j, "_",treatlevs[[j]][3],"v8_DEG.csv") )
      #write_csv(results(dds, contrast = c("treatment", treatlevs[[j]][2],treatlevs[[j]][3]))%>% 
      #            as.data.frame() %>% arrange(padj) %>% rownames_to_column('gene'),
      #          paste0('data/DESeq/', names(allcounts[i]), "_", j, "_",treatlevs[[j]][2],"v",treatlevs[[j]][3],"_DEG.csv") )
      write_csv(as.data.frame(counts(dds, normalized=T)) %>% rownames_to_column('gene'),
                paste0('data/normcounts/', names(allcounts[i]), "_normalized_counts.csv"))
      
      
      
      #print some stats along the way to make sure it's working
      print(results(dds, contrast = c("treatment", treatlevs[[j]][2], "8"))  %>% as.data.frame() %>% arrange(padj) %>% head(n=3))
      print(results(dds, contrast = c("treatment", treatlevs[[j]][2],"8"))  %>% as.data.frame() %>% filter(padj<0.05) %>% nrow())
      print(results(dds, contrast = c("treatment", treatlevs[[j]][3], "8"))  %>% as.data.frame() %>% arrange(padj) %>% head(n=3))
      print(results(dds, contrast = c("treatment", treatlevs[[j]][3],"8"))  %>% as.data.frame() %>% filter(padj<0.05) %>% nrow())
      
    }
  }
}
graph_single_gene('EPX', counts(dds, normalize=T), imeta)
graph_PCA(counts, imeta)
