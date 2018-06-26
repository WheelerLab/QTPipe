library(dplyr)
for (i in 1:22){

  sample_expressions <- as.data.frame(read.table(file = paste("/home/ryan/test_runs/juncfiles/Gene_expression_chr","1",sep=""), sep = '\t', header =T))
  reference_genotype <- as.data.frame(read.table(file = paste("/home/ryan/DATA/SNP_genotype_chr","1",".txt",sep=""), header=T))
  sampnames <- colnames(sample_expressions)
  refnames <- colnames(reference_genotype)
  filtnames <- base::intersect (sampnames, refnames)
  write.table(x = select(sample_expressions, "ID", filtnames), file = paste("/home/ryan/test_runs/juncfiles/filtered_Gene_expression_chr",i,sep=""), row.names = F, quote = F)
  write.table(x = select(reference_genotype, "rsid", filtnames), file = paste("/home/ryan/DATA/filtered_SNP_genotype_chr",i,".txt",sep=""), row.names = F, quote = F)  
  
}
  