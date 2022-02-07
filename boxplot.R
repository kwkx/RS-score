#URF-8 encoding, we used the codes below to draw boxplot for the internal rest set (GSE157103). If you want to draw 
#the box plots of other sets, please change the input file names to corresponding files.


setwd("D:\\新冠文章\\代码上传\\boxplot")#please set path to "boxplot" file

library(plyr)
library(ggplot2)
library(ggpubr)
scoreFile="157103_internal_test_RSscore.txt"    
cliFile="157103internal_testclinical.txt"            
trait="group"                    


score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(score), row.names(cli))
rt=cbind(score[sameSample,,drop=F], cli[sameSample,,drop=F])


bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,trait]))]
fix(rt)
rt[,c(trait, "group")]
colnames(rt1)=c("trait", "group")
df=as.data.frame(table(rt1))
#???y(df, .(group), transform, percent = Freq/sum(Freq) * 100)
#?ٷֱ?�
df=ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
df$group=factor(df$group, levels=c("Low", "High"))

�c(trait, "Score")]
colnames(rt2)=c("trait", "Score")
type=levels(factor(rt2[,"trait"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#???????t=
ggboxplot(rt2, x="trait", y="Score", fill="trait",
		          xlab=trait,
		          ylab="RS score",
		          legend.title=trait,
		          palette bioCol
		          )+  geom_jitter(width =0.2)+ 
	    stat_compare_means(comparisons=my_comparisons,method= "t.test")
pdf(file="boxplot.157103internal_test_boxplot.pdf4,height=4.5)
print(boxplot)
dev.off()
?geom_poi