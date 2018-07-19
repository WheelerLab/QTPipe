#install.packages("tidyr")
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))

parser <- ArgumentParser()
parser$add_argument("-a", "--annotation", help="file path of the annotation file")
parser$add_argument("-q", "--quantdir", help="file path of directory containing salmon's quantification files")
parser$add_argument("-o", "--outputdir", help="file path of the snp genotype file")
args <- parser$parse_args()

if ( ! exists(args$outputdir) ) {
  args$outputdir <- args$quantdir
}
#this chunk can be made more efficient
#also keep track of directories - these will be inputs for optparse
files <- list.files(path = args$quantdir , pattern = "[0123456789]+\\.sorted\\.genes\\.quant\\.sf", recursive = T)
sample_temp <- as.data.frame(read.table(file = paste(args$quantdir,files[1], sep=""), sep ='\t', header =T))
TPM <- select(sample_temp, "Name")
names <- c("gene_id", substring(files,1 ,regexpr("_", files)-1))
for (j in files){
  sample_temp <- as.data.frame(read.table(file = paste(args$quantdir,j, sep=""), sep='\t', header=T))
  TPM <- cbind(TPM, select(sample_temp, "TPM"))
}

colnames(TPM)<-names
print("TPM dataframe generated")

total_annotation <- as.data.frame(read.table(file = args$annotation , sep = '\t', header =F, skip = 5))
colnames(total_annotation) <- c("chr", "source", "feature-type", "start", "stop", "score", "strand", "phase", "metadata")
total_annotation<-tidyr::separate(total_annotation, metadata, c("gene_id", "transcript_id"), sep = ";", extra = "drop")
total_annotation$gene_id <- gsub("gene_id\\s(.*)", "\\1", total_annotation$gene_id)
total_annotation$transcript_id <- gsub("transcript_id\\s*(.*)", "\\1", total_annotation$transcript_id)
print("total annotation data frame generated")

for (i in 1:22){
  chr_annotation <- filter(total_annotation, chr == paste("chr", i, sep=""))
  write.table(x =select(semi_join(chr_annotation, TPM, by ="gene_id"), "gene_id", "chr", "start", "stop"), file = paste(args$outputdir, "/location_sal_chr",i,sep=""), row.names = F, quote = F)#location file
  write.table(x =semi_join(TPM, chr_annotation, by = "gene_id"),file = paste(args$outputdir, "/expression_sal_chr",i,sep=""), row.names = F, quote = F)
  print(paste("chr",i,"processed"))
}

