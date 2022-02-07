#utf-8 encoding
setwd("D:\\新冠文章\\代码上传\\ssGSEA")  #please change this path to the "ssGSEA" file.
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
####first, we perform ssGSEA to find out immune cell content in SHF and LHF groups, namelyu Figure 2A.
expFile="merge.txt"           
gmtFile="immunecontent.gmt"          
clusterFile="Cluster.txt"       
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA_immunecontent.result.txt",sep="\t",quote=F,col.names=F)

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

data=melt(scoreCluster, id.vars=c("cluster"))
colnames(data)=c("cluster", "Immune", "Fraction")

bioCol=c("#0066FF","#FF9900")
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="cluster", 
     ylab="Immune infiltration",
     xlab="",
     legend.title="cluster",
     palette=bioCol)
p=p+rotate_x_text(50)
pdf(file="immune_content_boxplot.pdf", width=8, height=6.5)                          #????ͼƬ?ļ?
p+stat_compare_means(aes(group=cluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()


####Then, we perform ssGSEA to find out immune function difference in SHF and LHF groups, namely Supplementary Figure 1.

expFile="merge.txt"           
gmtFile="immunefunction.gmt"          
clusterFile="Cluster.txt"       
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA_immunefunction.result.txt",sep="\t",quote=F,col.names=F)

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

data=melt(scoreCluster, id.vars=c("cluster"))
colnames(data)=c("cluster", "Immune", "Fraction")

bioCol=c("#0066FF","#FF9900")
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="cluster", 
            ylab="Immune infiltration",
            xlab="",
            legend.title="cluster",
            palette=bioCol)
p=p+rotate_x_text(50)
pdf(file="immune_function_boxplot.pdf", width=8, height=6.5)                          #????ͼƬ?ļ?
p+stat_compare_means(aes(group=cluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

