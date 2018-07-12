library(dplyr)
#argument 1 should be the fastq list with the names dict
#argument 2 should be the the SNP genotype txt file
arguments <- commandArgs(trailingOnly = T) 
sample_list <- read.table(arguments[1],sep='\t',header=F) #read in the samples
population_list <- read.table(arguments[2],header=F,nrows=1) #read in the population pool
population_list <- unlist(population_list[1,]) #convert population pool into vector
filtered_list <- filter(sample_list,V1 %in% population_list) #remove any samples that are not in the population pool of our snps
write.table(x = filtered_list, file = arguments[1], row.names = F, quote = F, col.names = F, sep = '\t') #write the new list to the old file
