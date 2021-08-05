#PCA主成分分析
library(ggfortify)
# 互换行和列，再dim一下
df=as.data.frame(t(b))
# 不要view df，列太多，软件会卡住；
dim(df)
dim(nor)

nor[1:6,1:6]
df[1:6,1:6]

df$group=group_list 
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')


a <- as.vector(c(x$study1,x$study2))
b <- as.vector(x$fisher)
length(a)
length(b)
a <- unique(a)
length(a)


b[which(b %in% a == FALSE)]




library(sva)
library(limma)
rt=read.table("merge.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
batchType=c(rep(1,50),rep(2,10))
modType=c(rep("normal",25),rep("tumor",25),rep("normal",5),rep("tumor",5))
mod = model.matrix(~as.factor(modType))
outTab=ComBat(data, batchType, mod, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab),outTab)
write.table(outTab,file="normalize.txt",sep="\t",quote=F,col.names=F)