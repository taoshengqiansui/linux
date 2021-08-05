#ENTREZID SYMBOL ENSEMBLE 的相互转换
suppressMessages(library(org.Hs.eg.db))# 载入包
keytypes(org.Hs.eg.db) #查看支持对选项
#rt=read.table("gene.txt",sep="\t",check.names=F,header=T)
rt <- pname$ENSEMBL
#rt<-as.data.frame(rt[,1])
class(rt)

#keys是自己的基因，columns是输出的类型，keytype是输入的类型
#gene_list<-select(org.Hs.eg.db, keys=as.character(rt$`rt[, 1]`), columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")
gene_list<-select(org.Hs.eg.db, keys=as.character(symbol), columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")
gene_list<-select(org.Hs.eg.db, keys=rt, columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")
gene_list <- gene_list[c(-3,-5,-60,-69,-93),]
gene_list[1:4,1:3]
write.table(gene_list,file="a.txt",sep="\t",quote=F,row.names=F)



#基因功能注释
##ID转换
library("org.Hs.eg.db")
rt <- pname$SYMBOL
rt <- na.omit(rt)
rt=read.table("symbol.txt",sep="\t",check.names=F,header=T)
genes=as.vector(gene_list[,1])
entrezIDs <- mget(rt, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
dim(out)
out<-out[(out$entrezID!='NA'),] #删除NA值
dim(out)

write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)

##2GO分析
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
rt=read.table("id.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]

gene=rt$entrezID

#GO富集分析
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)

#柱状图
tiff(file="GO.tiff",width = 26,height = 20,units ="cm",compression="lzw",bg="white",res=600)
barplot(kk, drop = TRUE, showCategory =10) #只展示一种
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#气泡图
tiff(file="dotplot.tiff",width = 26,height = 20,units ="cm",compression="lzw",bg="white",res=600)
dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#热图
tiff(file="heatplot.tiff",width = 40,height = 20,units ="cm",compression="lzw",bg="white",res=600)
heatplot(kk,showCategory =20, foldChange=cor)
dev.off()

#画柱状图
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  }
  else
  {
    return(x)
  }
}

data<-read.table("GO.txt",header = T,sep = "\t")
#替换名字
data$ONTOLOGY<-gsub("BP", "biological_process", data$ONTOLOGY)
data$ONTOLOGY<-gsub("CC", "cellular_component", data$ONTOLOGY)
data$ONTOLOGY<-gsub("MF", "molecular_function", data$ONTOLOGY)
class(data)
colnames(data)[3]<-c('GO_term')

data <- subset(data,Count>3)  #数目很多时才做,保证每个GO的基因数>3
#将GO_TERM的名字变短，shorten_names是上述定义的函数
data$GO_term=(sapply(levels(data$GO_term)[as.numeric(data$GO_term)],shorten_names))
data<-data[order(data[,1]),] #排序
data$GO_term<- as.character(data$GO_term)  #先转换成字符串
data$GO_term<-factor(data$GO_term,levels = c(data$GO_term)) #再强制加入因子

COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
a<-ggplot(data=data, aes(x=GO_term,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) + coord_flip() +  
  xlab("GO term") + ylab("Num of Genes") +
  scale_fill_manual(values = COLS)+ theme_bw()
ggsave(a, file="go_all.pdf", width=9.03, height=5.74)



#KEGG分析
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

rt=read.table("id.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]

gene=rt$entrezID


#kegg富集分析
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)

#柱状图
tiff(file="KEGG.tiff",width = 20,height = 12,units ="cm",compression="lzw",bg="white",res=600)
barplot(kk, drop = TRUE, showCategory = 20)
dev.off()

#气泡图
tiff(file="dotplot.tiff",width = 20,height = 12,units ="cm",compression="lzw",bg="white",res=600)
dotplot(kk, showCategory = 20)
dev.off()

#热图
tiff(file="heatplot.tiff",width = 25,height = 15,units ="cm",compression="lzw",bg="white",res=600)
heatplot(kk,showCategory =20, foldChange=cor)
dev.off()


#通路图
library("pathview")
keggxls=read.table("KEGG.txt",sep="\t",header=T)
for(i in keggxls$ID){
  pv.out <- pathview(gene.data = cor, pathway.id = i, species = "hsa", out.suffix = "pathview")
}
