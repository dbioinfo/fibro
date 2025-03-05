library(tidyverse)
library(data.table)
#list files
fnames <- list.files('data/counts/', recursive=T)

#read in to list of dataframes
all_abundances <- lapply(paste0('data/counts/',fnames), read_tsv)
names(all_abundances) <- sapply(str_split(basename(fnames), pattern='_'), '[', 1)

#calculate counts for each, format names better
counts <- imap(all_abundances, ~ .x %>%
  mutate(result = round((Coverage*(End-Start))/150,0)) %>%  
  rename('Transcript'=`Gene ID`) %>% 
  mutate(Transcript = sub("^([^.]+\\.[^.]+)\\..*$", "\\1", Transcript)) %>%
  group_by(Transcript) %>% #sometimes there's double values, sum them
  summarise(result = sum(result)) %>% 
  select(Transcript, !!.y := result) %>% #rename result to sample name
  ungroup()  
) 

# Get transcripts with at least 1 non-0 value in non-humans (human toga gtf unavailable)
unique_transcripts <- counts[36:length(counts)] %>%  
    map(~ filter(.x, .[[2]] != 0) %>% pull(Transcript)) %>%  
    reduce(union)
print(length(unique_transcripts))

#filter counts before merging
counts <- counts %>% 
    map(~ filter(.x, Transcript %in% unique_transcripts))
counts <- lapply(counts, as.data.table)

#merge dataframes (too much memory used here, have over 200Gb avail, still not enough, so each species)
hidx <- grepl(pattern='^HS', x=names(counts))
sqidx <- grepl(pattern='^OK|^FIT', x=names(counts))
lbidx <- grepl(pattern='^PG', x=names(counts))
fbidx <- grepl(pattern='^RA', x=names(counts))
ridx <- grepl(pattern='^RN', x=names(counts))
bidx <- grepl(pattern='^OK22', x=names(counts))

merged_counts <- reduce(counts[hidx], full_join, by='Transcript')
write_csv(merged_counts, 'data/human_full_count_matrix.csv')
print('human done')

merged_counts <- reduce(counts[sqidx & !bidx], full_join, by='Transcript')
write_csv(merged_counts, 'data/squirrel_full_count_matrix.csv')
print('squirrel done')

merged_counts <- reduce(counts[lbidx], full_join, by='Transcript')
write_csv(merged_counts, 'data/lbbat_full_count_matrix.csv')
print('lil brown bat done')

merged_counts <- reduce(counts[fbidx], full_join, by='Transcript')
write_csv(merged_counts, 'data/fruitbat_full_count_matrix.csv')
print('fruit bat done')

merged_counts <- reduce(counts[ridx], full_join, by='Transcript')
write_csv(merged_counts, 'data/rat_full_count_matrix.csv')
print('rat done')
