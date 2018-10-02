library(dplyr)
library(argparse)


parser <- ArgumentParser()
parser$add_argument("-s", "--samplelist", help="file path of the sample list")
parser$add_argument("-sg", "--snpgenotype", help="file path of the snp genotype file")
parser$add_argument("-sl", "--snplocation", help="file path of the snp location file")
parser$add_argument("-o", "--outputdir", help="output directory")
parser$add_argument("-t", "--tag", help="tag for this round of filtering")
args <- parser$parse_args()

SNP_genotype<-as.data.frame(read.table(file = args$snpgenotype, sep='', header = T))
print("loaded SNP data")
refnames <- colnames(SNP_genotype)
samplelist<-as.data.frame(read.table(file = args$samplelist, sep='\t', header = F))
if (ncol(samplelist) == 2){
sampnames<-make.names(samplelist[,2])
} else {
  sampnames<-make.names(samplelist[,1])
}
sampsubset<-base::intersect(refnames, sampnames)
samp_genotype<-select(SNP_genotype, "id", sampsubset)
samp_genotype<-samp_genotype[rowMeans(samp_genotype[-1]) != -1 
                           & rowMeans(samp_genotype[-1]) != 0 
                           & rowMeans(samp_genotype[-1]) != 1
                           & rowMeans(samp_genotype[-1]) != 2, ]
print("finished subsetting")
sampsubset<-gsub('.', '-', sampsubset, fixed = T)
SNP_location<-as.data.frame(read.table(file = args$snplocation, sep='', header = T))
write.table(x = semi_join(samp_genotype, SNP_location, by = c("id" = "snp")), file = paste(args$outputdir,"/",args$tag, "SNPgenotype.txt", sep = ''), row.names = F, quote = F)
write.table(x = semi_join(SNP_location, samp_genotype, by = c("snp" = "id")), file = paste(args$outputdir,"/",args$tag, "SNPlocation.txt", sep = ''), row.names = F, quote = F)
write.table(x=sampsubset, file = paste(args$outputdir,"/","Genotyped_subset.txt", sep = ''), row.names = F, quote = F, col.names = F)
