library(dplyr)
for (i in 1:22){

  sample_expressions <- as.data.frame(read.table(file = paste("/home/ryan/test_runs/expression_fiels/expression_sal_chr",i,sep=""), sep = '\t', header =T))
  reference_genotype <- as.data.frame(read.table(file = paste("/home/ryan/DATA/filtered_SNP_genotype_chr",i,".txt",sep=""), header=T))
  sampnames <- colnames(sample_expressions)
  refnames <- colnames(reference_genotype)
  filtnames <- base::intersect (sampnames, refnames)
  write.table(x = select(sample_expressions, "ID", filtnames), file = paste("/home/ryan/test_runs/expression_files/filtered_Gene_expression_chr",i,sep=""), row.names = F, quote = F)
}
