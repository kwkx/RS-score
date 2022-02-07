#UTF-8 encoding
setwd("D:\\新冠文章\\代码上传\\model")#please change the path to the "model" file

###############the codes for establishing the model(RS score) are shown below.
library("glmnet")
library("survival")
rt=read.table("training set.txt",header=T,sep="\t",row.names=1,check.names=F)     #读取文件
fix(rt)
x=as.matrix(rt[,c(2:ncol(rt))])
y=data.matrix(rt$label)
fix(y)
y=as.matrix(y)
data.class(y)

print(fit)
plot(fit, xvar="lambda", label=TRUE)
set.seed(112)
cvfit=cv.glmnet(x,y,type.measure="class",family="binomial" )
plot(cvfit)###we can get the CV plot here (Supplementary Figure 4A)
cvfit$lambda.min####lamda value:0.03111349

fit <- glmnet(x, y, family = "binomial", maxit = 1000)
plot(fit, xvar = "lambda", label = F) ###Supplementary Figure 4B

#####
fit <- glmnet(x, y, maxit = 1000,alpha=1,family = "binomial",lambda=0.03111349)#train the Lasso model(RS score)

l.coef1<-coef(fit,exact = T)
l.coef1=as.matrix(l.coef1)
write.table(l.coef1,file="coef.csv",sep=",",row.names=T,quote=F)
save(fit, file = "lassomodel.RData") 


###############We first test our RS score in the internal test set using below codes#########################
gene_list=read.table("gene_list.txt", header=T, sep="\t", check.names=F)
fix(gene_list)
pair_list=read.table("pair_list.txt", header=T, sep="\t", check.names=F)
fix(pair_list)
exp=read.table("testexp_beforepair.csv", header=T, sep=",", check.names=F,row.names = 1)
fix(exp)
out=data.frame()
for(i in 1:nrow(gene_list))
{
  num=which(row.names(exp)[]==gene_list[i,1])
  out=rbind(out,exp[num,])
  print(i)
  
}
#To pair the genes, we used codes below
tcgaPair=data.frame()
rt=out
fix(rt)
sampleNum=ncol(rt)
for(i in 1:(nrow(rt)-1)){
  for(j in (i+1):nrow(rt)){
    pair=ifelse(rt[i,]>rt[j,], 1, 0)
    rownames(pair)=paste0(rownames(rt)[i],"|",rownames(rt)[j])
    tcgaPair=rbind(tcgaPair, pair)
  }
}

tcgaOut=rbind(ID=colnames(tcgaPair), tcgaPair)
tcgaOut=tcgaOut[which(rownames(tcgaOut)[]%in%pair_list[,1]),]
fix(tcgaOut)
nrow(tcgaOut)


tcgaOut=t(tcgaOut)
tcgaOut=as.data.frame(tcgaOut)
for (i in 1:ncol(tcgaOut)) {
  tcgaOut[,i]=as.numeric(tcgaOut[,i])
}
tcgaOut=as.matrix(tcgaOut)
riskScoreTest=predict(fit, newx = tcgaOut[,])      #"fit" is the RS score model.
out=cbind(row.names(tcgaOut),riskScoreTest)
fix(out)
label=read.table("testlabel.csv", header=T, sep=",", check.names=F,row.names = 1)#"1" means LHF group, "0" means SHF group.
fix(label)
out=cbind(out,label)
colnames(out)=c("ID","Score","label")

write.table(out,
            file="internal_test_set_score.csv",
            sep=",",
            quote=F,
            row.names=F)
###This file can be used to creat ROC curve on "http://www.rocplot.org/site/index"
###1.choosing "generate ROC plot using uploaded data"
###2.upload the "internal_test_set_score.csv" file
###3."select marker" chooses "Score";"Select status" chooses "label"; "Response equals" chooses "1";"Separator" chooses "Comma"
###4."Run"




###############We then test our RS score in the GSE177477 cohort using below codes#########################
gene_list=read.table("gene_list.txt", header=T, sep="\t", check.names=F)
fix(gene_list)
pair_list=read.table("pair_list.txt", header=T, sep="\t", check.names=F)
fix(pair_list)
exp=read.table("GSE177477_beforepair.csv", header=T, sep=",", check.names=F,row.names = 1)
fix(exp)
out=data.frame()
for(i in 1:nrow(gene_list))
{
  num=which(row.names(exp)[]==gene_list[i,1])
  out=rbind(out,exp[num,])
  print(i)
  
}
#To pair the genes, we used codes below
tcgaPair=data.frame()
rt=out
fix(rt)
sampleNum=ncol(rt)
for(i in 1:(nrow(rt)-1)){
  for(j in (i+1):nrow(rt)){
    pair=ifelse(rt[i,]>rt[j,], 1, 0)
    rownames(pair)=paste0(rownames(rt)[i],"|",rownames(rt)[j])
    tcgaPair=rbind(tcgaPair, pair)
  }
}

tcgaOut=rbind(ID=colnames(tcgaPair), tcgaPair)
tcgaOut=tcgaOut[which(rownames(tcgaOut)[]%in%pair_list[,1]),]
fix(tcgaOut)
nrow(tcgaOut)


tcgaOut=t(tcgaOut)
tcgaOut=as.data.frame(tcgaOut)
for (i in 1:ncol(tcgaOut)) {
  tcgaOut[,i]=as.numeric(tcgaOut[,i])
}
tcgaOut=as.matrix(tcgaOut)
riskScoreTest=predict(fit, newx = tcgaOut[,])      #"fit" is the RS score model.
out=cbind(row.names(tcgaOut),riskScoreTest)
label=read.table("GSE177477label.csv", header=T, sep=",", check.names=F,row.names = 1)#"1" means Symptomatic group, "0" means Asymptomatic group.
finalname=intersect(row.names(out),row.names(label))
out=out[finalname,]
label=label[finalname,]
out=cbind(out,label)
colnames(out)=c("ID","Score","label")
write.table(out,
            file="GSE177477_score.csv",
            sep=",",
            quote=F,
            row.names=F)
###This file can be used to creat ROC curve on "http://www.rocplot.org/site/index"
###1.choosing "generate ROC plot using uploaded data"
###2.upload the "GSE177477_score.csv" file
###3."select marker" chooses "Score";"Select status" chooses "label"; "Response equals" chooses "1";"Separator" chooses "Comma"
###4."Run"



###############We then test our RS score in the GSE172114 cohort using below codes#########################
gene_list=read.table("gene_list.txt", header=T, sep="\t", check.names=F)
fix(gene_list)
pair_list=read.table("pair_list.txt", header=T, sep="\t", check.names=F)
fix(pair_list)
exp=read.table("GSE172114_beforepair.csv", header=T, sep=",", check.names=F,row.names = 1)
fix(exp)
out=data.frame()
for(i in 1:nrow(gene_list))
{
  num=which(row.names(exp)[]==gene_list[i,1])
  out=rbind(out,exp[num,])
  print(i)
  
}
#To pair the genes, we used codes below
tcgaPair=data.frame()
rt=out
fix(rt)
sampleNum=ncol(rt)
for(i in 1:(nrow(rt)-1)){
  for(j in (i+1):nrow(rt)){
    pair=ifelse(rt[i,]>rt[j,], 1, 0)
    rownames(pair)=paste0(rownames(rt)[i],"|",rownames(rt)[j])
    tcgaPair=rbind(tcgaPair, pair)
  }
}

tcgaOut=rbind(ID=colnames(tcgaPair), tcgaPair)
tcgaOut=tcgaOut[which(rownames(tcgaOut)[]%in%pair_list[,1]),]
fix(tcgaOut)
nrow(tcgaOut)


tcgaOut=t(tcgaOut)
tcgaOut=as.data.frame(tcgaOut)
for (i in 1:ncol(tcgaOut)) {
  tcgaOut[,i]=as.numeric(tcgaOut[,i])
}
tcgaOut=as.matrix(tcgaOut)
riskScoreTest=predict(fit, newx = tcgaOut[,])      #"fit" is the RS score model.
out=cbind(row.names(tcgaOut),riskScoreTest)
label=read.table("GSE172114label.csv", header=T, sep=",", check.names=F,row.names = 1)#"1" means Critical group, "0" means Non-critical group.
finalname=intersect(row.names(out),row.names(label))
out=out[finalname,]
label=label[finalname,]
out=cbind(out,label)
colnames(out)=c("ID","Score","label")
write.table(out,
            file="GSE172114_score.csv",
            sep=",",
            quote=F,
            row.names=F)
###This file can be used to creat ROC curve on "http://www.rocplot.org/site/index"
###1.choosing "generate ROC plot using uploaded data"
###2.upload the "GSE172114_score.csv" file
###3."select marker" chooses "Score";"Select status" chooses "label"; "Response equals" chooses "1";"Separator" chooses "Comma"
###4."Run"




###############We then test our RS score in the GSE155454 cohort using below codes#########################
gene_list=read.table("gene_list.txt", header=T, sep="\t", check.names=F)
fix(gene_list)
pair_list=read.table("pair_list.txt", header=T, sep="\t", check.names=F)
fix(pair_list)
exp=read.table("GSE155454_beforepair.csv", header=T, sep=",", check.names=F,row.names = 1)
fix(exp)
out=data.frame()
for(i in 1:nrow(gene_list))
{
  num=which(row.names(exp)[]==gene_list[i,1])
  out=rbind(out,exp[num,])
  print(i)
  
}
#To pair the genes, we used codes below
tcgaPair=data.frame()
rt=out
fix(rt)
sampleNum=ncol(rt)
for(i in 1:(nrow(rt)-1)){
  for(j in (i+1):nrow(rt)){
    pair=ifelse(rt[i,]>rt[j,], 1, 0)
    rownames(pair)=paste0(rownames(rt)[i],"|",rownames(rt)[j])
    tcgaPair=rbind(tcgaPair, pair)
  }
}

tcgaOut=rbind(ID=colnames(tcgaPair), tcgaPair)
tcgaOut=tcgaOut[which(rownames(tcgaOut)[]%in%pair_list[,1]),]
fix(tcgaOut)
nrow(tcgaOut)


tcgaOut=t(tcgaOut)
tcgaOut=as.data.frame(tcgaOut)
for (i in 1:ncol(tcgaOut)) {
  tcgaOut[,i]=as.numeric(tcgaOut[,i])
}
tcgaOut=as.matrix(tcgaOut)
riskScoreTest=predict(fit, newx = tcgaOut[,])      #"fit" is the RS score model.
out=cbind(row.names(tcgaOut),riskScoreTest)
label=read.table("GSE155454label.csv", header=T, sep=",", check.names=F,row.names = 1)#"1" means severity 2 group, "0" means severity 0+1 group.
finalname=intersect(row.names(out),row.names(label))
out=out[finalname,]
label=label[finalname,]
out=cbind(out,label)
colnames(out)=c("ID","Score","label")
write.table(out,
            file="GSE155454_score.csv",
            sep=",",
            quote=F,
            row.names=F)
###This file can be used to creat ROC curve on "http://www.rocplot.org/site/index"
###1.choosing "generate ROC plot using uploaded data"
###2.upload the "GSE155454_score.csv" file
###3."select marker" chooses "Score";"Select status" chooses "label"; "Response equals" chooses "1";"Separator" chooses "Comma"
###4."Run"



###############We then test our RS score in the GSE152418 cohort using below codes#########################
gene_list=read.table("gene_list.txt", header=T, sep="\t", check.names=F)
fix(gene_list)
pair_list=read.table("pair_list.txt", header=T, sep="\t", check.names=F)
fix(pair_list)
exp=read.table("GSE152418_beforepair.csv", header=T, sep=",", check.names=F,row.names = 1)
fix(exp)
out=data.frame()
for(i in 1:nrow(gene_list))
{
  num=which(row.names(exp)[]==gene_list[i,1])
  out=rbind(out,exp[num,])
  print(i)
  
}
#To pair the genes, we used codes below
tcgaPair=data.frame()
rt=out
fix(rt)
sampleNum=ncol(rt)
for(i in 1:(nrow(rt)-1)){
  for(j in (i+1):nrow(rt)){
    pair=ifelse(rt[i,]>rt[j,], 1, 0)
    rownames(pair)=paste0(rownames(rt)[i],"|",rownames(rt)[j])
    tcgaPair=rbind(tcgaPair, pair)
  }
}

tcgaOut=rbind(ID=colnames(tcgaPair), tcgaPair)
tcgaOut=tcgaOut[which(rownames(tcgaOut)[]%in%pair_list[,1]),]
fix(tcgaOut)
nrow(tcgaOut)


tcgaOut=t(tcgaOut)
tcgaOut=as.data.frame(tcgaOut)
for (i in 1:ncol(tcgaOut)) {
  tcgaOut[,i]=as.numeric(tcgaOut[,i])
}
tcgaOut=as.matrix(tcgaOut)
riskScoreTest=predict(fit, newx = tcgaOut[,])      #"fit" is the RS score model.
out=cbind(row.names(tcgaOut),riskScoreTest)
label=read.table("GSE152418label.csv", header=T, sep=",", check.names=F,row.names = 1)#"1" means ICU group, "0" means Moderate + Severe group.
finalname=intersect(row.names(out),row.names(label))
out=out[finalname,]
label=label[finalname,]
out=cbind(out,label)
colnames(out)=c("ID","Score","label")
write.table(out,
            file="GSE152418_score.csv",
            sep=",",
            quote=F,
            row.names=F)
###This file can be used to creat ROC curve on "http://www.rocplot.org/site/index"
###1.choosing "generate ROC plot using uploaded data"
###2.upload the "GSE152418_score.csv" file
###3."select marker" chooses "Score";"Select status" chooses "label"; "Response equals" chooses "1";"Separator" chooses "Comma"
###4."Run"



