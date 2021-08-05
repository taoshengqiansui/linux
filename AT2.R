library(metaRNASeq)
study2 <- read.table(file = "./GSE152586.matrix",header = T,row.names = 1)
study4 <- read.table(file = "./GSE155518.matrix",header = T,row.names = 1)
study4 <- study4[,c(4:6,1:3)]
study6 <- read.table(file = "./GSE160435.matrix",header = T,row.names = 1)
study6 <- study6[,c(3,4,6,9,10,1,2,5,7,8)]
study4 <- study4[rownames(study2),]
study6 <- study6[rownames(study2),]

c1 <- c(rep("Normal",times = 5),rep("Infected",times = 16))
c2 <- c(rep("Normal",times = 3),rep("Infected",times = 3))
c3 <- c(rep("Normal",times = 3),rep("Infected",times = 3))
c4 <- c(rep("Normal",times = 3),rep("Infected",times = 3))
c6 <- c(rep("Normal",times = 5),rep("Infected",times = 5))
c7 <- c(rep("Normal",times = 3),rep("Infected",times = 3))
colData1 <- data.frame("study"=c(rep("study1",times = 21)),"condition"=c1,
                       row.names = c(paste("Normal",c(1:5),sep = ""),paste("Infected",c(1:16),sep = "")))
colData2 <- data.frame("study"=c("study2","study2",
                                 "study2","study2","study2","study2"),
                       "condition"=c2,row.names = c("Normal1","Normal2","Normal3",
                                                    "Infected1","Infected2","Infected3"))
colData3 <- data.frame("study"=c(rep("study3",times = 6)),"condition"=c3,
                       row.names = c(paste("Normal",c(1:3),sep = ""),
                                     paste("Infected",c(1:3),sep = "")))
colData4 <- data.frame("study"=c(rep("study4",times = 6)),"condition"=c4,
                       row.names = c(paste("Normal",c(1:3),sep = ""),paste("Infected",c(1:3),sep = "")))
colData6 <- data.frame("study"=c("study6","study6","study6","study6","study6","study6",
                                 "study6","study6","study6","study6"),
                       "condition"=c6,row.names = c(paste("Normal",c(1:5),sep = ""),
                                                    paste("Infected",c(1:5),sep = "")))
colData7 <- data.frame("study"=c(rep("study7",times = 6)),
                       "condition"=c7,row.names = c(paste("Normal",c(1:3),sep = ""),
                                                    paste("Infected",c(1:3),sep = "")))
if (requireNamespace("DESeq2", quietly = TRUE)) {
  dds2 <- DESeq2::DESeqDataSetFromMatrix(countData = study2,
                                         colData = colData2,design = ~ condition)
  res2 <- DESeq2::results(DESeq2::DESeq(dds2))
  
  dds4 <- DESeq2::DESeqDataSetFromMatrix(countData = study4,
                                         colData = colData4,design = ~ condition)
  res4 <- DESeq2::results(DESeq2::DESeq(dds4))
  
  dds6 <- DESeq2::DESeqDataSetFromMatrix(countData = study6,
                                         colData = colData6,design = ~ condition)
  res6 <- DESeq2::results(DESeq2::DESeq(dds6))}


  rawpval <- list("pval1"=res2[["pvalue"]],"pval2"=res4[["pvalue"]],"pval3"=res6[["pvalue"]])
  FC <- list("FC1"=res2[["log2FoldChange"]],"FC2"=res4[["log2FoldChange"]],"FC3"=res6[["log2FoldChange"]])

##  获取每个数据的Fold change和原始P值并合并为Fold change列表和P值列表
  adjpval <- list("adjpval1"=res2[["padj"]],"adjpval2"=res4[["padj"]],"adjpval3"=res6[["padj"]])
## 获取每个数据的矫正后的P值并合并为列表
studies <- c("study1","study2","study3")
DE <- mapply(adjpval, FUN=function(x) ifelse(x <= 0.05, 1, 0))
colnames(DE)=paste("DE",studies,sep=".")
#studies为前面定义的一个含有两个字符串的向量
#三.检验数据是否符合均匀分布，即检验合并P值的方法是否可行
par(mfrow = c(2,2))
hist(rawpval[[1]], breaks=100, col="grey", main="Study 1", xlab="Raw p-values")
hist(rawpval[[2]], breaks=100, col="grey", main="Study 2", xlab="Raw p-values")
hist(rawpval[[3]], breaks=100, col="grey", main="Study 3", xlab="Raw p-values")

##结果不是很好，因为rawpval中有很多的NA(DESeq2的分析结果)，不确定
filtered <- lapply(adjpval, FUN=function(pval) which(is.na(pval)))
rawpval[[1]][filtered[[1]]]=NA
rawpval[[2]][filtered[[2]]]=NA
rawpval[[3]][filtered[[3]]]=NA
par(mfrow = c(2,2))
hist(rawpval[[1]], breaks=100, col="grey", main="Study 1",xlab="Raw p-values")
hist(rawpval[[2]], breaks=100, col="grey", main="Study 2",xlab="Raw p-values")
hist(rawpval[[3]], breaks=100, col="grey", main="Study 3", xlab="Raw p-values")

##结果好了很多
#四.meta合并
##1.fisher合并
metacomb <- fishercomb(rawpval, BHth = 0.05)
hist(metacomb$rawpval, breaks=100, col="grey", main="Fisher method",
     xlab = "Raw p-values (meta-analysis)")
##2.逆正太合并
#invnormcomb <- invnorm(rawpval,nrep=c(8,8), BHth = 0.05)
#hist(invnormcomb$rawpval, breaks=100, col="grey",
#main="Inverse normal method",
#xlab = "Raw p-values (meta-analysis)")
DEresults <- data.frame(DE,
                        "DE.fishercomb"=ifelse(metacomb$adjpval<=0.05,1,0))
#,"DE.invnorm"=ifelse(invnormcomb$adjpval<=0.05,1,0))
#DE已经有两列,再加两列,DEresults只显示了是否是差异表达基因
head(DEresults)
signsFC <- mapply(FC, FUN=function(x) sign(x))
##sign()函数结果:正为1,0为0,负为-1.
sumsigns <- apply(signsFC,1,sum)
##margin=1/2,1代表行,2代表列.
##如果signsFC中有表达不一致的基因,则结果为0.
##所以这一步是筛选出表达一致的(都是上调或下调)的基因.
commonsgnFC <- ifelse(abs(sumsigns)==dim(signsFC)[2], sign(sumsigns),0)
##除了表达一致的基因为(2变为1)1,其他都为0
unionDE <- unique(c(metacomb$DEindices))
#invnormcomb$DEindices))
##筛选出两种合并方法中的所有差异基因所对应的索引,包括只在一种方法中显示的差异基因
##分别length():fishcomb$DEindices,invnormcomb$DEindices,unionDE,结果为1452,1390,1460
FC.selecDE <- data.frame(DEresults[unionDE,],do.call(cbind,FC)[unionDE,],
                         signFC=commonsgnFC[unionDE])
FC.selecDE$DE <- ifelse(abs(FC.selecDE$signFC)==1,TRUE,FALSE)
FC.selecDE <- FC.selecDE[!is.na(FC.selecDE$DE),]
##去除DE=NA所在的行
##do.call传递命令,命令为cbind,对象为FC,FC按列合并,也就是将两个数据的合并成了两列
##4+2+1+1
keepDE <- FC.selecDE[which(abs(FC.selecDE$signFC)==1),]
##两套数据所得到的表达一致的基因保留下来
conflictDE <- FC.selecDE[which(FC.selecDE$signFC == 0),]
##表达矛盾的基因
dim(FC.selecDE)
dim(keepDE)
dim(conflictDE)
head(keepDE)
table(conflictDE$DE)
#韦恩图
fishcomb_de <- rownames(keepDE)[which(keepDE[,"DE.fishercomb"]==1)]
#invnorm_de <- rownames(keepDE)[which(keepDE[,"DE.invnorm"]==1)]
indstudy_de <- list(rownames(keepDE)[which(keepDE[,"DE.study1"]==1)],
                    rownames(keepDE)[which(keepDE[,"DE.study2"]==1)],
                    rownames(keepDE)[which(keepDE[,"DE.study3"]==1)])

length(fishcomb_de)
#length(invnorm_de)
length(indstudy_de[[1]])
length(indstudy_de[[2]])
length(indstudy_de[[3]])

IDD.IRR(fishcomb_de,indstudy_de)
#IDD.IRR(invnorm_de ,indstudy_de)
if (require("VennDiagram", quietly = TRUE)) {
  venn.plot<-venn.diagram(x = list(study1=which(keepDE[,"DE.study1"]==1),
                                   study2=which(keepDE[,"DE.study2"]==1),
                                   study3=which(keepDE[,"DE.study3"]==1),
                                   fisher=which(keepDE[,"DE.fishercomb"]==1)),
                          filename = NULL, col = "black",
                          fill = c("blue","purple","green","yellow"),
                          margin=0.05, alpha = 0.6)
  jpeg("at2_jpeg.jpg");
  grid.draw(venn.plot);dev.off();}
sessionInfo()


##所有DEG
###ID转换
x = list(study1=which(keepDE[,"DE.study1"]==1),
         study2=which(keepDE[,"DE.study2"]==1),
         study3=which(keepDE[,"DE.study3"]==1),
         study4=which(keepDE[,"DE.study4"]==1),
         study6=which(keepDE[,"DE.study6"]==1),
         fisher=which(keepDE[,"DE.fishercomb"]==1))
all.ensemblID <- rownames(study1[as.integer(rownames(keepDE)),])
fisher <- x$fisher
fisherDE <- keepDE[,c(-1:-6)]

FC_all <- data.frame(fisherDE$FC1,fisherDE$FC2,fisherDE$FC3,fisherDE$FC4,fisherDE$FC6)

fisherDE$FC <- rowMeans(FC_all)
fisherDE <- fisherDE[,c(6,8)]
rownames(fisherDE) <- all.ensemblID
library("org.Hs.eg.db")
meta_symbol<-select(org.Hs.eg.db, keys=rownames(fisherDE), columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")
rownames(meta_symbol) <- meta_symbol$ENSEMBL# 找出重复ID
meta_symbol <- meta_symbol[c(-9,-577,-384,-18,-314,-685,-418,-353),]
meta_symbol <- cbind(meta_symbol,fisherDE)
length(which(is.na(meta_symbol$SYMBOL)))#90个NA
meta.symbol <- meta_symbol[which(is.na(meta_symbol$SYMBOL)==FALSE),]

library(clusterProfiler)
library(ggplot2)
ego <- enrichGO(gene = meta.symbol$ENTREZID, OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", ont = "ALL",qvalueCutoff=0.05)

barplot(ego,showCategory = 8,split ="ONTOLOGY")+facet_grid(ONTOLOGY~.,scales = "free")
ega <- as.data.frame(ego)
View(ega)

kk <- enrichKEGG(gene = meta.symbol$ENTREZID,organism = "hsa",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
kk <- setReadable(kk,org.Hs.eg.db,keyType = "ENTREZID")
allgene.fc <- meta.symbol$FC
names(allgene.fc) <- meta.symbol$SYMBOL
allgene.fc <- sort(allgene.fc,decreasing = T)

barplot(kk,showCategory = 10)
dotplot(kk,showCategory = 10)
cnetplot(kk,foldChange = allgene.fc,showCategory = 10,node_label = "all",circular = TRUE,colorEdge =TRUE, color_gene = "red",cex_label_gene = 0.3)
#不可用emapplot(kk,showCategory = 6,pie = "count")






#FC排列
meta_symbol$signFC <- ifelse(meta_symbol$signFC==1,"UP","DOWN")
rm.na <- meta_symbol[which(is.na(meta_symbol$SYMBOL)==FALSE),]
all.rmna.FC <- head(as.data.frame(rm.na[order(abs(rm.na$FC),decreasing = T),]),50)
all.order.FC <- head(as.data.frame(meta_symbol[order(abs(meta_symbol$FC),decreasing = T),]),50)
a <- meta_symbol[which(meta_symbol$ENSEMBL %in% meta$name==TRUE),]

##上下调的基因symbol
library(clusterProfiler)
library(ggplot2)
up.ID <- rm.na[which(rm.na$signFC=="UP"),]$ENTREZID
down.ID <- rm.na[which(rm.na$signFC=="DOWN"),]$ENTREZID
kk.up <- enrichKEGG(gene = up.ID,organism = "hsa",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
barplot(kk.up,showCategory = 10)
dotplot(kk.up,showCategory = 10)
cnetplot(kk.up,foldChange = allgene.fc,showCategory = 10,node_label = "gene",circular = TRUE,colorEdge =TRUE)
kk.down<- enrichKEGG(gene = up.ID,organism = "hsa",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
barplot(kk.down,showCategory = 10)
dotplot(kk.down,showCategory = 10)
cnetplot(kk.down,foldChange = allgene.fc,showCategory = 10,node_label = "gene",circular = TRUE,colorEdge =TRUE)

order.FC <- head(as.data.frame(a[order(abs(a$FC),decreasing = T),]),20)
save(meta_symbol,meta.symbol,all.order.FC,pname,all.ensemblID,all.rmna.FC,order.FC,file = "./meta.Rdata")




##新发现DEG
x = list(study1=which(keepDE[,"DE.study1"]==1),
         study2=which(keepDE[,"DE.study2"]==1),
         study3=which(keepDE[,"DE.study3"]==1),
         study4=which(keepDE[,"DE.study4"]==1),
         study6=which(keepDE[,"DE.study6"]==1),
         fisher=which(keepDE[,"DE.fishercomb"]==1))
a <- sort(as.vector(c(x$study1,x$study2,x$study3,x$study4,x$study6)))
b <- as.vector(x$fisher)
length(a)
length(b)
a <- unique(a)
length(a)
list <- b[which(b %in% a == FALSE)]
#list指的是keepDE的索引而不是行名
keepde <- keepDE[list,]
list2 <- as.integer(rownames(keepDE[list,]))
#list2指的是study1中的索引，而不是行名
name <- rownames(study1[list2,])


p <- as.data.frame(metacomb$rawpval)
p <- p[list2,]
adjp <- as.data.frame(metacomb$adjpval)
adjp <- adjp[list2,]
meta <- data.frame(name,p,adjp,keepde$signFC)
meta$DE <- ifelse(meta$keepde.signFC==1,"UP","DOWN")



FC2 <- data.frame(keepde$FC1,keepde$FC2,keepde$FC3,keepde$FC4,keepde$FC6)
meta$FC <- rowMeans(FC2)

pmeta <- cbind(meta,pname)
pmeta <- pmeta[which(is.na(pmeta$SYMBOL)==F),]
genelist <- pmeta[,c("SYMBOL","FC")]
fc <- as.vector(genelist$FC)
names(fc) <- genelist$SYMBOL

gene_list <- gene_list[c(-3,-5,-60,-69,-93),]
glist <- gene_list[is.na(gene_list$ENTREZID)==FALSE,]

