library(dplyr)

#Subset SNP samples
SNP_genotype<-as.data.frame(read.table(file = '/home/ryan/Data/SNPdata/GEU_SNPs.txt', sep='\t', header = T))
print("loaded SNP data")
refnames <- colnames(SNP_genotype)
samplelist<-as.data.frame(read.table(file = '/home/ryan/sample_list.txt', sep='\t', header = F))
sampnames<-samplelist[,1]
sampsubset<-base::intersect(refnames, sampnames)
samp_genotype<-select(SNP_genotype, "rsid", sampsubset)
print("finished subsetting")
for (i in 1:22){
  SNP_location<-as.data.frame(read.table(file = paste('/home/ryan/Data/SNPdata/SNP_Location_chr',i,'.txt', sep =''), sep='\t', header = T))
  write.table(x = semi_join(samp_genotype, SNP_location), file = paste('/home/ryan/Data/SNPdata/SNP_genotype_chr',i,'.txt', sep=''))
  print(paste("chr",i,"processed"))
}
