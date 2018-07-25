library(dplyr)

#Subset SNP samples
SNP_genotype<-as.data.frame(read.table(file = '/home/ryan/Data/SNPdata/GEU_SNPs.txt', sep='\t', header = T))
print("loaded SNP data")
refnames <- colnames(SNP_genotype)
samplelist<-as.data.frame(read.table(file = '/home/ryan/sample_list.txt', sep='\t', header = F))
sampnames<-samplelist[,1]
sampsubset<-base::intersect(refnames, sampnames)
samp_genotype<-select(SNP_genotype, "rsid", sampsubset)
samp_genotype<-samp_genotype[rowMeans(samp_genotype[-1]) != 0 
                           & rowMeans(samp_genotype[-1]) != 1
                           & rowMeans(samp_genotype[-1]) != 2, ]
print("finished subsetting")
for (i in 1:22){
  SNP_location<-as.data.frame(read.table(file = paste('/home/ryan/Data/SNPdata/SNP_Location_chr',i,'.txt', sep =''), sep='\t', header = T))
  write.table(x = semi_join(samp_genotype, SNP_location), file = paste('/home/ryan/Data/SNPdata/SNP_genotype_chr',i,'.txt', sep=''), row.names = F, quote = F)
  write.table(x = semi_join(SNP_location, samp_genotype), file = paste('/home/ryan/Data/SNPdata/SNP_Location_chr',i,'.txt', sep =''), row.names = F, quote = F)
  print(paste("chr",i,"processed"))
}
