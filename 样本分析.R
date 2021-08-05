#样本相关性分析
study1 <- read.table(file = "./GSE155518.matrix",header = T,row.names = 1)
cor1 <- round(cor(study1),digits = 2)
#cor() 函数可以指定，method="kendall"/"spearman"
#默认是皮尔森（线性相关），斯皮尔曼（等级相关），肯德尔（适用于离散型/分类型变量）
library(pheatmap)
pheatmap(cor1)


#样本聚类
##计算距离矩阵
dist1 <- dist(t(study1))

##聚类
hc1 <- hclust(dist1)
plot(hc1)


#主成分分析
library(PCAtools)
pca(study1)




library(clusterProfiler)
MF <- enrichGO(OrgDb = "org.Hs.eg.db",gene = rownames(a),pvalueCutoff = 0.05,readable = TRUE)



#GO
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(topGO)
library(pathview)
library(Rgraphviz)
enid <- rownames(a)
id <- mapIds(x = org.Hs.eg.db,keys = enid,keytype = "ENSEMBL",column="ENTREZID")
di <- na.omit(id)
length(di)
BP <- enrichGO(gene = di,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = T)



#UpSetR
library(UpSetR)
venn <- list(study_1=x$study1,study_2=x$study2,study_3=x$study3,study_4=x$study4,study_5=x$study6,meta_analysis=x$fisher)
color <- c("#EA4335","#FBBC05","#34A853","#4285f4","#4DAF4A","#FF7F00")
UpSetR::upset(fromList(venn),nsets = 6,sets = "A",point.size = 3,line.size = 0,
              order.by = "degree",matrix.color = "black",
              mainbar.y.label = "Intersection Size",
              sets.bar.color = color,
              queries = list(list(query = intersects,
                                  param =list("A","B","C"),
                                  color = "orange",active = T)))

#nset=6显示所有这六个数据
#matrix.color指的是交集的颜色，也就是下方的原点
#main.bar.color 表示柱子的颜色;sets.bar.color左下角柱子颜色
#mainbar.y.labe表示y轴标签
#order.by = "freq"按照上图柱子长度排列，order.by = "degree，按照下图原点数排列


upset(fromList(venn), order.by = "degree",nsets = 6,
      point.size = 3,matrix.color = "orange",mainbar.y.label = "Intersection Size",
    main.bar.color = "yellow",sets.bar.color = "blue"
    list(list(query=intersects,params = list("meta_analysis"),color = "orange"))
    )
upset(fromList(venn), order.by = "degree",nsets = 6,
      point.size = 3,matrix.color = "orange",mainbar.y.label = "Intersection Size",
      main.bar.color = "black",sets.bar.color = "blue",
      queries = list(list(query=intersects,params = list("meta_analysis"),
                          color = "red")))
