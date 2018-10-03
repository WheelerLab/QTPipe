library(dplyr)
library(argparse)

parser$add_argument("-s", "--samplelist", help="file path of the sample list")
parser$add_argument("-sg", "--snpgenotype", help="file path of the snp genotype file")
parser$add_argument("-o", "--outputdir", help="output directory", type=)

sample_list <- read.table(args$samplelist,sep='\t',header=F) #read in the samples
population_list <- read.table(args$snpgenotype,header=F,nrows=1) #read in the population pool
population_list <- unlist(population_list[1,]) #convert population pool into vector
if (ncol(samplelist) == 2){
  filtered_list <- filter(sample_list,V2 %in% population_list) #remove any samples that are not in the population pool of our snps
} else {
  filtered_list <- filter(sample_list,V1 %in% population_list) #remove any samples that are not in the population pool of our snps
}
write.table(x = filtered_list, file = paste(args$outputdir, "/genotyped_samples.txt"), row.names = F, quote = F, col.names = F, sep = '\t') #write the new list to the old file
