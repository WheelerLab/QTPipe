"%&%" = function(a,b) paste(a,b,sep="")
library(peer)
library(tibble)
library(dplyr)
source('~/scripts/rntransform.R')
source('~/scripts/ztransform.R')
library(preprocessCore)
library(argparse)
library(methods)
parser <- ArgumentParser()
parser$add_argument("-e", "--expression", help="file path of the gene expression file")
parser$add_argument("-o", "--outputdir", help="output directory", type="character")
parser$add_argument("-t", "--tag", help="file tag for this run of samples")
parser$add_argument("--cov", help="file of covariates")
args <- parser$parse_args()
expr = read.table(file=args$expression,header=T,check.names = F)

expmat = as.matrix(t(expr[,-1])) #need genes in cols, ids in rows
qn.expmat <- as.matrix(normalize.quantiles(expmat)) ##quantile normalize
qn<-as.matrix(t(qn.expmat))
rownames(qn) <- expr[,1]
colnames(qn) <- colnames(expr)[-1]
write.table(file = "~/mets_analysis/normalized_expression/" %&% args$tag %&% ".txt", x = qn, quote = F, sep = '\t')

#rn.qn.expmat <- apply(qn.expmat,1,"rntransform") ##rank transform to normality & transposes##
#rownames(rn.qn.expmat) <- expr[,1]
#colnames(rn.qn.expmat) <- colnames(expr)[-1]
#write.table(file = "~/mets_analysis/normalized_expression/rank_normalized_splicing_chr" %&% args$tag %&% ".txt.gz", x = rn.qn.expmat, quote = F, sep = '\t')
#covs = as.matrix(read.table(args$cov,header=T, row.names = 1))
#covs<-as.matrix(as.numeric(covs))
  ## i is number of peer factors to calculate, recommend 25% of sample size, but no more than 100
for(i in c(10,15,20)){
    model = PEER()
    PEER_setPhenoMean(model, qn.expmat)
    #PEER_setCovariates(model, as.matrix(covs))
    PEER_setNk(model,i)
    PEER_update(model)
    factors = PEER_getX(model)
    weights = PEER_getW(model)
    precision = PEER_getAlpha(model)
    residuals = PEER_getResiduals(model)
    print(dim(residuals))
    cat(args$expression, "Nk =", i)
    print(PEER_plotModel(model))
    adjexp <- t(residuals)
    print(dim(adjexp))
    
    rownames(adjexp) <- expr[,1]
    colnames(adjexp) <- colnames(expr)[-1]
    adjexp <- as.data.frame(adjexp) %>% rownames_to_column('gene_id')
    png(filename =args$outputdir %&% "/" %&% args$tag %&% "_" %&% i %&% '_peer.png')
    factors = PEER_getX(model)
    rownames(factors) <- colnames(expr)[-1]
    write.table(factors,file=args$outputdir %&% "/" %&% args$tag %&% "_splicing_expression" %&% i %&% ".PEER.factors" %&% ".txt", quote=F)
    PEER_plotModel(model)
    dev.off()
    write.table(adjexp, file=args$outputdir %&% "/" %&% args$tag %&% "_splicing_expression" %&% i %&% '_peer_factor_adjusted.txt' , quote=F, row.names=F, sep='\t')
}


