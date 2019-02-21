library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
FDRcounts<-as.data.frame(read.table(file="Z:/mets_analysis/meqtl/lauren_ALL/fdr_counts.txt", header = F, stringsAsFactors = F))
#FixedQNcounts<-as.data.frame(read.table(file="Z:/mets_analysis/meqtl/FDR_fixed_QN_counts.txt", header = F))
#FDRcounts<-filter(FDRcounts, V1 != 0)


FDRcounts["CHR"]<-NA
FDRcounts["Cov"]<-NA
FDRcounts["PF"]<-NA
FDRcounts["normalization"]<-NA
FDRcounts["ID"]<-NA
#highest<-rep(NA,22)
#FixedQNcounts["CHR"]<-NA
#FixedQNcounts["Cov"]<-NA
#FixedQNcounts["PF"]<-NA
#FixedQNcounts["normalization"]<-NA
#FixedQNcounts["ID"]<-NA
#FDRcounts<-FDRcounts[c(1:330),]
for (i in c(1:88)){
  string<-FDRcounts[i,2]
  FDRcounts[i,3]<-str_extract(string, "chr[0-9]+")
  if (grepl("Nk_10",string) || grepl("pf10",string)){
    FDRcounts[i,5]<-10
  } else if (grepl("Nk_15",string) || grepl("pf15",string)){
    FDRcounts[i,5]<-15
  } else if (grepl("Nk_20",string) || grepl("pf20",string)){
    FDRcounts[i,5]<-20
  } else if (grepl("Nk_30",string) || grepl("pf30",string)){
    FDRcounts[i,5]<-30
  } else {
    FDRcounts[i,5]<-0
  }
    
  if (grepl("pc10",string)){
    FDRcounts[i,4]<-10
  } else if (grepl("pcs_3",string)){
    FDRcounts[i,4]<-3
  } else {
    FDRcounts[i,4]<-0
  }
  FDRcounts[i,6]<-if_else(grepl("rank_normalized" ,string),true = "QN.RN",false = "QN")
  FDRcounts[i,7]<-paste("PC",FDRcounts[i,4],"PF",FDRcounts[i,5],FDRcounts[i,6],sep=".")
}
#FDRcounts<-filter(FDRcounts, PF != 0)
#FDRcounts<-rbind(FDRcounts,FixedQNcounts)
#write.table(x=FDRcounts,file = "Z:/mets_analysis/meqtl_gene_counts_pop_correction.txt", col.names = T, sep = '\t', row.names = F, quote = F)
#dat<-separate(FDRcounts,V2,sep="_", into = c(NA,NA,V1,V2,PF))
#new<-chr[which.max(chr$V1),]
#for (i in c(1:22)) {
#  chr<-filter(FDRcounts,CHR==paste("chr",i,sep=""))
#  new<-rbind(new,chr[which.max(chr$V1),])
  
#}
#?geom_bar
#p<-ggplot(data = parameter_counts, aes(x=parameters, y = count))
#p+geom_col()
#test<-subset(counts, Cov == "3" & CHR == "chr1")
#which.max(chr$V1)
new2<-unique(FDRcounts$ID)
parameter_counts<-data.frame(count=rep(NA,4),parameters=rep(NA,4))
for (i in c(1:4)){

  parameterized<-filter(FDRcounts, ID==new2[i])
  if (nrow(parameterized) != 22) {
    print("fuck")
  }
  parameter_counts[i,1]<-sum(parameterized$V1)
  parameter_counts[i,2]<-new2[i]
}
parameter_counts["PF"]<-NA
parameter_counts["PC"]<-NA
parameter_counts["normalization"]<-NA
for (i in c(1:4)){
  string<-parameter_counts[i,2]
  parameter_counts[i,3]<-str_extract(string, "PF\\.[0-9]+")
  parameter_counts[i,4]<-str_extract(string, "PC\\.[0-9]+")
  parameter_counts[i,5]<-if_else(grepl("RN" ,string),true = "QN.RN",false = "QN")
}
#?gsub()
#quantile_counts<-filter(parameter_counts,normalization!="QN.RN")
#parameter_counts$PF<-gsub(x=parameter_counts$PF,pattern="PF.",fixed=T,replacement = "")
p2<-ggplot(data = parameter_counts, aes(x = PF,y=count))# + geom_point(aes(colour=factor(PF)))
p2 + geom_point(aes(colour=factor(PC))) + geom_line(aes(group=factor(PC), colour=factor(PC))) +ylab("cisQTLs at FDR < 0.05") + ggtitle("PF vs QTL count for Lauren's ALL population") #+ ylim(0,470000)
#p3<-ggplot(data = quantile_counts, aes(x = PC,y=count))
#p3 + geom_point(aes(colour=factor(PF))) + geom_line(aes(group=factor(PF), colour=factor(PF)))
