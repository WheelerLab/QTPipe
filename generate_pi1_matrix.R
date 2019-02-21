library(data.table)
library(dplyr)
library(ggplot2)
library(qvalue)
library(argparse)
"%&%" = function(a,b) paste(a,b,sep="")

parser <- ArgumentParser()
#parser$add_argument("--discovery", help="path to the discovery file")
parser$add_argument("--replication", help="path to the replication file")
parser$add_argument("--out", help="file you would like to output as")
args <- parser$parse_args()

#List METS and MESA pops
discovery_cohorts<-c("METS_WG_FDR0.05_PC0_PF0.txt.gz","METS_WG_FDR0.05_PC0_PF10.txt.gz","METS_WG_FDR0.05_PC0_PF15.txt.gz","METS_WG_FDR0.05_PC0_PF20.txt.gz","METS_WG_FDR0.05_PC10_PF0.txt.gz","METS_WG_FDR0.05_PC10_PF10.txt.gz","METS_WG_FDR0.05_PC10_PF15.txt.gz","METS_WG_FDR0.05_PC10_PF20.txt.gz","METS_WG_FDR0.05_PC3_PF0.txt.gz","METS_WG_FDR0.05_PC3_PF10.txt.gz","METS_WG_FDR0.05_PC3_PF15.txt.gz","METS_WG_FDR0.05_PC3_PF20.txt.gz")
#replication_cohorts<-as.matrix(read.table(file="~/mets_analysis/meqtl/replication_pops/pop_list"))
discovery_dir<-"~/mets_analysis/meqtl/all_cis/"
#replication_dir<-"~/mets_analysis/meqtl/replication_pops/"

#initialize an empty pi1 Matrix
pi1_matrix<-matrix(NA, nrow = length(discovery_cohorts), ncol = 1)
rownames(pi1_matrix)<-discovery_cohorts
colnames(pi1_matrix)<-args$replication

#read in replication files
print("reading in replication file")
replication_file<-fread(file=args$replication, header = T,stringsAsFactors = F, colClasses = "character", sep = "\t", showProgress = T)
colnames(replication_file)<-c("snps","gene","statistic","pvalue","FDR","beta")

print("finished reading in replication file")
for (i in c(1:length(discovery_cohorts))){
    #read in discovery files
    print("reading in discovery file")
    discovery_file<-fread("zcat " %&% discovery_dir %&% discovery_cohorts[i], stringsAsFactors = F, colClasses = "character", sep = "\t", header = T, showProgress = T)
    colnames(discovery_file)<-c("snps","gene","statistic","pvalue","FDR","beta", "FDR_WG")
    discovery_file$gene<-gsub("\\..*", "", discovery_file$gene)
    print("finished reading in discovery file")
    row<-discovery_cohorts[i]
    #select the  snp-gene pairs that replicate
    replicated<-inner_join(discovery_file, replication_file, by=c("snps","gene"))
    str(replicated)
    #grab statistics for pi1 calculation
    over<-dim(replicated)
    #teat<-junda[,9]==0#lapply(junda[,9], function(x) x[is.infinite(x)])
    #print(teat)
    replicated_pvals<-as.numeric(replicated$pvalue.y)
    png("~/pvals_ALL_rep.png")
    hist(replicated_pvals)
    dev.off()
    str(replicated_pvals)
    #print("STR of bad0")
    #bad0<-replicated_pvals[replicated_pvals==0]
    #str(bad0)
    #print("str of badinf")
    #badinf<-lapply(replicated_pvals, function(x) x[is.infinite(x)])
    #str(badinf)
    #calculate pi1
    qobj_replicated<-qvalue(p = replicated_pvals)
    pi1<- 1 - qobj_replicated$pi0
    #truncate pi1
    pi12<- signif(pi1,4)
    col<-args$replication
    #add to matrix
    pi1_matrix[row,col] <- pi12
    cat("finished combination " %&% i %&% " " %&% args$replication)
}

write.table(pi1_matrix, file=args$out,quote=F, row.names = T, sep = '\t')