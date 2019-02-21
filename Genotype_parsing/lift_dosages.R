## script for converting dosage files to build hg38
##uses hg19_hg38_rsid_dict_wg.txt, which drops alternative chromosomes
#dosage files are hardcoded in as input for translation
library(data.table)
library(dplyr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-o", "--outputdir", help="where to output")
#parser$add_argument("-d", "--dict", help="dictionary .txt file to be read in")
args <- parser$parse_args()

"%&%" = function(a,b) paste(a,b,sep="")

dosage_dir<-"~/mets_analysis/genotypes_postimputation/sanger_imputation/Sanger_dosages/"
print("reading in dictionary file")
dict<-fread(dosage_dir %&% "hg19_hg38_rsid_dict_wg.txt", header = T, sep = "\t", showProgress = T)
#dict<-fread(dosage_dir %&% "test.dict.txt", header = T, sep = "\t", showProgress = T)
print("finished reading in dict")
for (i in c(1:22)){
  print("reading in chr " %&% i %&% " dosages")
  dosage<-fread("zcat " %&% dosage_dir %&% "chr" %&% i %&% ".maf0.01.info0.8.dosage.txt.gz")
  #dosage<-fread(dosage_dir %&% "test.dosage.txt", header = T)
  print("done")
  print("original number of snps: " %&% dim(dosage)[1])
  names<-colnames(dosage)
  names[2]<-"hg38.cpos"
  names[3]<-"hg38.pos"
  print("performing inner join")
  joined.df<-inner_join(dosage, dict, by = c("snp_ID" = "hg19.cpos"))
  joined.df<-select(joined.df, names)
  print("done")
  print("final number of snps: " %&% dim(joined.df)[1])
  names[2]<-"snp_ID"
  names[3]<-"pos"
  colnames(joined.df)<-names
  print("writing to file")
  fwrite(joined.df, file = args$outputdir %&% "/hg38.chr" %&% i %&% ".maf0.01.info0.8.dosage.txt", col.names = T, row.names = F, sep = " ", quote = F)
}
print("Done. Have a nice day :)")
