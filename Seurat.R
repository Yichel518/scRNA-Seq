setwd("~/scRNA-Seq/new_10X/New_0.8/")
## Single cell QC
metadata <-ovary@meta.data
# 为元数据添加细胞ID
metadata$cells <- rownames(metadata)
# 细胞计数 (Cell counts)
# 每个细胞的UMI计数 (UMI counts per cell)
# 每个细胞检测到的基因 (Genes detected per cell)
# 检测到的UMI数对比基因数 (UMIs vs. genes detected)
# 线粒体计数比率 (Mitochondrial counts ratio)
# 复杂度 (Novelty)
# 重命名列
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ovary"))] <- "ovary"
#metadata$sample[which(str_detect(metadata$cells, "^testis"))] <- "testis"
metadata$sample <- metadata$seq_folder
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
# 可视化每个细胞的UMI/转录本数目
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# 通过频数图可视化每个细胞检测出的基因数分布
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# 通过箱线图可视化每个细胞检测到的基因的分布
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
# 可视化检测到的基因数和UMI数之间的关系，并且观察是否存在大量低数目的基因数/UMI数的细胞
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=percent.mito)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
# 可视化每个细胞检测到的线粒体基因表达分布
metadata %>% 
  ggplot(aes(color=sample, x=percent.mito, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
# 通过可视化每一个UMI检测到的基因数来可视化基因表达的整体复杂性
ovary$log10GenesPerUMI <- log10(ovary$nFeature_RNA) / log10(ovary$nCount_RNA)
log10GenesPerUMI <- ovary$log10GenesPerUMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
# 1.Setup the Seurat Object（创建Seurat对象）
library(dplyr)
library(Seurat)
rm(list = ls())
ovar.data <- Seurat::Read10X(data.dir = "~/scRNA-Seq/zebrafish_OVAR/")
##########TESTDATA#############################
te = load10X(system.file('extdata','toyData',package='SoupX'))
mDat = read.table(system.file('extdata','toyData','metaData.tsv',package='SoupX'),sep='\t')
te = setClusters(te,mDat$res.1)
te = autoEstCont(te)
out = adjustCounts(te)
#############################
tsne_ovary=ovary@reductions$tsne@
cluster_matrix=as.data.frame(ovary$RNA_snn_res.0.3)
metadata=cbind(cluster_matrix,tsne_ovary)
#write.table(metadata,"~/scRNA-Seq/ovary_outs/metadata.tsv",col.names=T,row.names=T,quote=F,sep="\t")

sc = load10X("~/scRNA-Seq/ovary_outs/")
metadata=as.data.frame(ovary@meta.data)
sc=setClusters(sc, setNames(metadata$seurat_clusters, 
                                         rownames(metadata)))
sc = autoEstCont(sc)
out = adjustCounts(sc)
ovary_data=out
# # 第一点：基因在多少细胞表达 
#fivenum(apply(ovar.data,1,function(x) sum(x>0) ))

#boxplot(apply(ovar.data,1,function(x) sum(x>0) ))

# 第二点：细胞中有多少表达的基因
#fivenum(apply(ovar.data,2,function(x) sum(x>0) ))
#hist(apply(ovar.data,2,function(x) sum(x>0) ))
#dim(ovar.data)


# Initialize the Seurat object with the raw (non-normalized data).
ovary <- CreateSeuratObject(counts = ovar.data, project = "ovary", min.cells = 3, 
                            min.features =200) # gene:13714,cells：2700
# 比较稀疏矩阵和普通矩阵使用的内存大小
dense.size <- object.size(x = as.matrix(x = ovar.data))
dense.size
sparse.size <- object.size(x = ovar.data)
sparse.size
dense.size / sparse.size

# 2.Standard pre-processing workflow（数据预处理流程）
ovary[["percent.mito"]] <- PercentageFeatureSet(object = ovary, pattern = "^mt-") #线粒体基因以MT开头的选择出来
ovary.mt=ovary[["percent.mito"]]
## nFeature_RNA就是每个细胞检测到的基因数目；nCount_RNA就是这些基因数目一共测到的count数目,也就是总的UMI数目；percent.mt（使用'PercentageFeatureSet`函数计算线粒体QC指标，该函数计算源自一组特征的计数百分比）
# Show QC metrics for the first 5 cells
head(x = ovary@meta.data, 5)

# Visualize QC metrics as a violin plot 用小提琴图可视化QC特征
VlnPlot(object = ovary, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3) 

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# FeatureScatter通常用于可视化特征 - 特征关系

plot1 <- FeatureScatter(object = ovary, feature1 = "nCount_RNA", feature2 = "percent.mito") 
plot2 <- FeatureScatter(object = ovary, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
CombinePlots(plots = list(plot1,plot2))
# 选择gene数目小于5000 以及大于200 以及 线粒体数目小于5%的细胞 #可以基于violin图来选择你要卡的阈值
# 过滤完之后还剩余 genes:13714,cells：2638
ovary <- subset(x = ovary, subset = nFeature_RNA > 200 & nFeature_RNA <5000& percent.mito < 20)
#（2）数据标准化
# 使用log转化，度量单位是10000
ovary <- NormalizeData(object = ovary, normalization.method = "LogNormalize", scale.factor = 1e4)

#（3）Identification of highly variable features (feature selection) 识别高度变异基因(特征选择的过程)、
ovary <- FindVariableFeatures(object = ovary,selection.method = 'dispersion', nfeatures = 2500)
# Identify the 10 most highly variable genes 取前十个变化最大的基因
top20 <- head(x = VariableFeatures(object = ovary), 20)
# plot variable features with and without labels 画出2500个高变异基因
plot1 <- VariableFeaturePlot(object = ovary)
LabelPoints(plot = plot1, points = top20, repel = TRUE)

#（4）Scaling 中心化
#all.genes <- rownames(x = ovary)
#ovary <- ScaleData(object = ovary, features = all.genes)
#- `ScaleData` function to remove unwanted sources of variation from a single-cell dataset, 比如mitochondrial contamination线粒体污染
ovary <- ScaleData(object = ovary, vars.to.regress = c("nCount_RNA", "percent.mito"))

# 3.Perform linear dimensional reduction线性降维 + Cluster the cells 对细胞进行聚类分析
##接下来，我们对中心化的数据执行PCA。 默认情况下，只有先前确定的变量要素用作输入，但如果要选择其他子集，则可以使用`features`参数定义。

#ovary <- SCTransform(ovary, vars.to.regress = "percent.mito", verbose = FALSE)


ovary <- RunPCA(object = ovary, features = VariableFeatures(object = ovary))

# Examine and visualize PCA results a few different ways                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
#print(x = ovary[['pca']], dims = 1:5, nfeatures = 5)
#VizDimLoadings(object = ovary, dims = 1:2, reduction = 'pca')
#DimPlot(object = ovary, reduction = 'pca') #每个细胞在PC1 和 PC2中的位置图
# Doheatmap 画热图
#DimHeatmap(object = ovary, dims = 1, cells = 500, balanced = TRUE) #可以选择所要呈现的细胞数目
#DimHeatmap(object = ovary, dims = 1:3, cells = 500, balanced = TRUE) 
# Determine the 'dimensionality' of the dataset 决定选取多少主成分
# 使用jackStraw方法 来估计多少主成分比较合适
# We identify 'significant' PCs as those who have a strong enrichment of low p-value features. 显著的主成分是强烈富集 低P值的特征
ovary <- JackStraw(object = ovary, num.replicate = 100)
ovary <- ScoreJackStraw(object = ovary, dims = 1:20)
# JackStrawplot
JackStrawPlot(object = ovary, dims = 1:20) #JackStrawPlot功能提供了一个可视化工具，用于比较每个PC的p值分布和均匀分布
# 碎石图
ElbowPlot(object = ovary) #基于每一个解释的方差百分比（“ElbowPlot”函数）对主成分进行排序

ovary <- FindNeighbors(object = ovary, dims = 1:16)
ovary <- FindClusters(object = ovary, resolution = c(seq(0.4,1.6,0.2)))
clustree(ovary@meta.data,prefix="RNA_snn_res.")
# Look at cluster IDs of the first 5 cells
head(x = Idents(object = ovary), 5)
# 4.Run non-linear dimensional reduction (UMAP/tSNE) 非线性降维： tSNE and UMAP
#reticulate::py_install(packages = "umap-learn") #安装UMAP
# UMAP
ovary=readRDS("~/scRNA-Seq/new_10X/ovary20200909.rds")
ovary <- FindClusters(object = ovary, resolution =0.8)
ovary <- RunUMAP(object = ovary, dims = 1:16)
#t=as.data.frame(ovary@active.ident)
#data=cbind(t,u)
#u=ovary@reductions$umap@cell.embeddings
#delet=data[data$`ovary@active.ident`==11,]%>%subset(.,UMAP_1<0)
#ovary@meta.data=ovary@meta.data[-which(rownames(ovary@meta.data)%in%rownames(delet)),]

DimPlot(object = ovary, reduction = 'tsne',label = TRUE)
DimPlot(object = ovary, reduction = 'umap',label = TRUE,
        cols=c("#88ada6","#D8BFD8","#800000","#ffc773","#a1afc9",
               "#4682B4","#808080","#556B2F","#D2B48C",
               "#5F9EA0","#DAA520","#BDB76B","#F4A460"),pt.size = 0.8)

##"#2E8B57",,"#A0522D","#C0C0C0","#FFB6C1","#4169E1","#008B8B","#FF6347"
# TSNE
ovary <- RunTSNE(object = ovary, dims = 1:16)
DimPlot(object = ovary, reduction = 'tsne',label = TRUE )
DimPlot(object = ovary, reduction = 'tsne',label = TRUE,
        cols=c("#2E8B57","#88ada6","#D8BFD8","#800000","#ffc773","#a1afc9",
               "#4682B4","#808080","#556B2F","darkgreen","#D2B48C",
               "#5F9EA0","#DAA520","#BDB76B","#F4A460","#A0522D",
               "#C0C0C0","#FFB6C1","#4169E1","#008B8B","#FF6347"),pt.size = 0.8)

#ovary <- RunICA(object = ovary)
#DimPlot(object = ovary, reduction = 'ica',label = TRUE )
# 您可以在此时保存对象，以便可以轻松地将其重新加载，而无需重新运行上面执行的计算密集型步骤，或者可以轻松地与协作者共享
# 5. Finding differentially expressed features 寻找差异表达的特征（聚类生物标志物）
# #- (1)找出所有cluster中的marker
# find all markers of cluster 1 找出cluster1 的所有marker
#cluster1.markers <- FindMarkers(object = ovar, ident.1 = 1, min.pct = 0.25)
#head(x = cluster1.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3 找到cluster 5 和cluster(0,3)比较的所有marker
#cluster5.markers <- FindMarkers(object = ovar, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(x = cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones 找所有cluster的marker
ovary_new.markers <- FindAllMarkers(object = ovary, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(ovary_new.markers,"~/scRNA-Seq/new_10X/New_0.8/ovary_new.markers20210204.txt",col.names=T,row.names=T,quote=F,sep="\t")
top50_ovar.new=ovary_new.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) #将ovar.markers传给group_by，group_by按cluster 排序，再将结果传给top_n，top_n按avg_logFC排序，显示每个类中的前两个
write.table(top50_ovar.new,"~/scRNA-Seq/new_10X/New_0.8/top50_ovar20210204.txt",col.names=T,row.names=T,quote=F,sep="\t")
top5_ovar.new=ovary_new.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) #将ovar.markers传给group_by，group_by按cluster 排序，再将结果传给top_n，top_n按avg_logFC排序，显示每个类中的前两个
write.table(top5_ovar.new,"~/scRNA-Seq/new_10X/New_0.8/top5_ovar20210204.txt",col.names=T,row.names=T,quote=F,sep="\t")
##################################################


#- (2)可视化这些marker在这些集群中的表达情况

#--1）用VlnPlot(显示跨集群的表达概率分布)
#--2）FeaturePlot（在tSNE或PCA图上可视化特征表达）这两种是比较常用的
#--3）还可以用`RidgePlot`，`CellScatter`和`DotPlot`这几种方法查看数据集情况

# VlnPlot
#ovary=readRDS("/home/tuyx/scRNA-Seq/new_10X/ovary20200909.rds")
DotPlot(object = ovary, features = c("dnd1","ddx4","sycp3","piwil1"))

DotPlot(object = ovary, features = c("piwil1","ddx4","sycp3","zp3.2","zar1"))
VlnPlot(object = OVARY, features = c("sycp3"))
DotPlot(object = ovary, features = c("zp3.2","zar1"))
DotPlot(object = ovary, features = c("eed","suz12a","phc2a"))
# FeaturePlot
FeaturePlot(object = ovary, features = c("zp3.2","zar1"))
FeaturePlot(object = ovary, features = c("bmi1a",
                                         "eed",
                                         "ezh2",
                                         "nanos2",
                                         "phc2a",
                                         "suz12a"
                                    ))
FeaturePlot(object = ovary, features = c("foxp3a","ddx4"
))
#RidgePlot
RidgePlot(object = ovar, features = c("bmi1a",
                                      "eed",
                                      "ezh2",
                                      "nanos2",
                                      "phc2a",
                                      "suz12a","ddx4","sycp3"))

#CellScatter
CellScatter(object = ovary, features = c("cyp17a1"),cell1 = "CATTACACGGAGTG", cell2 = "CATTACACCAACTG")

# DotPlot
DotPlot(object = ovary, features = c("cyp17a1","gsdf","cyp19a1a",
                                   "angptl4","fcer1g","grn1","spi1b",
                                   "fli1a","sox7","ccr9a","coro1a",
                                   "ddx4","piwil1","sycp3","col1a1b",
                                   "krt15","zar1","zp3b","cyp26a1"))
##oogonia7/4
# DotPlot(object =ovary, features = c("bmi1a",
#                                      "eed",
#                                      "ezh2",
#                                      "nanos2",
#                                      "phc2a",
#                                      "suz12a"
#                                     ))

##Epithelia 3 #6
DotPlot(object = ovary, features = c("ddx4","dnd1","nanos3","ca15b","h1m","rgs14a","slc16a3","wee2"))
## liver cell
DotPlot(object = ovary, features = c("fabp10a","tfa"))
## mast cell
DotPlot(object =ovary, features = c("tnfrsf18","ca2","sla2","rgs13","tnfrsf9b","cxcr4b"))
##Mast cells ok
DotPlot(object = ovary, features = c("fcer1g"))
DotPlot(object =ovary, features = c("fcer1g","cpa5","gata2","myd88"))
## oocyte
DotPlot(object = ovary, features = c("sycp3","sycp1","smc1b"))
VlnPlot(object = ovary, features = c("fcer1g"))
##Macrophages ok 

DotPlot(object = ovary, features = c("grn1","spi1b"))

VlnPlot(object = ovary, features = c("grn1","spi1b"))
VlnPlot(object = OVARY, features = c("grn1","spi1b"))
VlnPlot(object = Germ_cell, features = c("grn1","spi1b"))
##Granulosa 4 #3
DotPlot(object = ovary, features = c("cyp26a1","gsdf","bmp15","cyp19a1a","gata4","pgr","rspo1"))
DotPlot(object =ovary, features = c("cyp26a1","gsdf","angpt14","cyp19a1a","amh"))
##Vasculature ok #15
DotPlot(object = ovary, features = c("fli1a","sox7"))
VlnPlot(object = OVARY, features = c("fli1a","sox7"))
VlnPlot(object = ovary, features = c("fli1a","sox7"))
##Germ cells
###6
DotPlot(object = ovary, features = c(ovary_stage$I))
DotPlot(object = Germ_cell, features = c(ovary_stage$I))
###
DotPlot(object = ovary, features = c(ovary_stage$II))
DotPlot(object = ovary, features = ovary_stage$III)
DotPlot(object = ovary, features = ovary_stage$IV)
##meiosis
DotPlot(object = ovary, features = c("cdc25b","cdc25d","birc5b","aurka"))
##Sertoli cells 9 #8
DotPlot(object = ovary, features = c("cyp26a1"))
VlnPlot(object = OVARY, features = c("cyp26a1"))
VlnPlot(object =Germ_cell, features = c("cyp26a1"))
##Immuse cells
DotPlot(object = ovary, features = c("ccr9a","coro1a"))
VlnPlot(object = ovary, features = c("ccr9a","coro1a"))
##Connective tissue  ok  #2
DotPlot(object = ovary, features = c("col1a1b","krt15"))
VlnPlot(object = OVARY, features = c("col1a1b","krt15"))
##Theca cells 6
DotPlot(object = ovary, features = c("cyp17a1","star","bmp15"))
VlnPlot(object = OVARY, features = c("cyp17a1","star","bmp15"))
DotPlot(object = ovary, features = c("dnd","nanos3","ca15b","h1m","rgs14a","slc16a3","wee2","tdrd7a","dazl","gra","ddx4"))
DotPlot(object = ovary, features = c("dnd1","nanos3","ca15b","h1m","rgs14a","slc16a3","wee2"))

DotPlot(object = ovary, features = c("ddx4","sycp3","pou5f1"))

new.cluster.ids <- c("Oocyte I.o","Oocyte I.o","Granulosa cell.o","Lymphocyte-like.o",
                     "Lymphocyte-like.o", "Gonadal somatic cell.o","Oocyte I.o","Sertoli-like cell.o",
                     "Oogonia.o","Myeloblast-like.o","Gonadal somatic cell.o","Gonadal somatic cell.o",
                     "B cell-like.o","Granulosa cell.o","Erythroid cell.o","Endothelial cell.o",
                     "Oogonia.o","Regulatory T cell-like.o", "Hepatocyte.o","Myeloblast.o")
names(x=new.cluster.ids)=levels(x = ovary)
ovary <- RenameIdents(object =ovary, new.cluster.ids)
df1=as.data.frame(ovary@active.ident)
ovary@meta.data$celltype=df1$`ovary@active.ident`
saveRDS(ovary,"~/scRNA-Seq/new_10X/New_0.8/ovary20210204.rds")
####cell correlation
ovary<-CalculateBarcodeInflections(ovary)
SubsetByBarcodeInflections(ovary)
AverageExp.o<-AverageExpression(ovary,features=unique(top50_ovar.new$gene))
library(psych)
library(pheatmap)
coorda<-corr.test(AverageExp.o$RNA,AverageExp.o$RNA,method="spearman")
pheatmap(coorda$r,border_color = F,fontsize_col = 16,fontsize_row = 16)
######Biological processes#########################
cell <- unique(c("Oocyte I.o","Oocyte I.o","Granulosa cell.o","Lymphocyte-like.o",
                     "Lymphocyte-like.o", "Gonadal somatic cell.o","Oocyte I.o","Sertoli-like cell.o",
                     "Oogonia.o","Myeloblast-like.o","Gonadal somatic cell.o","Gonadal somatic cell.o",
                     "B cell-like.o","Granulosa cell.o","Erythroid cell.o","Endothelial cell.o",
                     "Oogonia.o","Regulatory T cell-like.o", "Hepatocyte.o","Myeloblast.o"))
setwd("~/scRNA-Seq/new_10X/New_0.8/BP_ovary")
for (i in 1:length(cell)){
  fileName=unlist(paste0(cell[i],".","txt"))
  tmp1=ovary_new.markers[ovary_new.markers$cluster==cell[i],]
  tmp2<- AnnotationDbi::select(org.Dr.eg.db,keys=tmp1$gene,columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
  BP.o<- enrichGO(gene=tmp2$SYMBOL,keyType = "SYMBOL",OrgDb= org.Dr.eg.db, ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.05)
  BP.o <- simplify(BP.o)
  BP.o_t <- as.data.frame(BP.o)
  write.table(BP.o_t,sprintf("%s",fileName),col.names=T,row.names=T,quote=F,sep="\t")
} 
setwd("~/scRNA-Seq/new_10X/New_0.8/KEGG_ovary")
for (i in 1:length(cell)){
  fileName=unlist(paste0(cell[i],".","txt"))
  tmp1=ovary_new.markers[ovary_new.markers$cluster==cell[i],]
  tmp2<- AnnotationDbi::select(org.Dr.eg.db,keys=tmp1$gene,columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
  path <- enrichKEGG(gene =tmp2$ENTREZID,
                              organism = 'dre',
                              pvalueCutoff = 0.05)
  path_table <- as.data.frame(path)
  write.table(path_table,sprintf("%s",fileName),col.names=T,row.names=T,quote=F,sep="\t")
} 


library(tidyverse)
library(magrittr)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(aplot)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
# 加载前面预处理好的scRNA-Seq数据，要用到表达矩阵、cluster及配色

# 每个基因在每个cluster里的平均值
avgData <- ovary@assays$RNA@data[top50_ovar.new$gene,] %>% 
  apply(1, function(x){
    tapply(x, ovary$celltype, mean) # ExpMean
  }) 
head(avgData)
phData <- MinMax(t(scale(avgData)), -2, 2) # z-score
#rownames(phData) <- 1:nrow(phData)
head(phData)
gene_features <- top50_ovar.new$gene
cluster_info <- sort(ovary$celltype)
col_o=data.frame(celltype=levels(ovary$celltype),color=c("#88ada6","#D8BFD8","#800000","#ffc773","#a1afc9",
                                                              "#4682B4","#808080","#556B2F","#D2B48C",
                                                              "#5F9EA0","#DAA520","#BDB76B","#F4A460"))
col=c("#88ada6","#D8BFD8","#800000","#ffc773","#a1afc9",
      "#4682B4","#808080","#556B2F","#D2B48C",
      "#5F9EA0","#DAA520","#BDB76B","#F4A460")

names(col) <- levels(cluster_info)

Heatmap(phData,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = FALSE,
        show_row_names = FALSE)

#为了增加聚类注释，我们需要用到HeatmapAnnotation函数，它对细胞的列进行注释，
#而rowAnnotation函数可以对行进行注释。这两个函数能够增加各种类型的注释，
#包括条形图，点图，折线图，箱线图，密度图等等，这些函数的特征是anno_xxx，
#例如anno_block就用来绘制区块图。

#top_anno <- HeatmapAnnotation(
#  cluster = anno_block(gp = gpar(fill = col), # 设置填充色
#                       labels = NULL, #levels(cluster_info), 
#                       labels_gp = gpar(cex = 0.8, col = "white"))) # 设置字体
top_anno = HeatmapAnnotation(df = data.frame(CellType = levels(ovary$celltype)),
                              col = list(CellType = col), show_legend = F)

select <- phData[rownames(phData)%in%c("fn1b","cxcl12a","cd81a","gsdf","myl9a","sox4a","fabp10a","tfa","leg1.1","csf3a",
                                   "cebpa","spi1b","lyve1a","mcamb","sox18","gnai1","tnfb","ccr9a","rac2",
                                   "grap2a", "nfkbiab","hbaa1","hbaa2","hbba2","jun","junba","ddx4","il7r","foxp3a",
                                   "piwil1","dnd1","zp3","nanos3","ddx4","sycp1","buc","fabp10a","serpina1","fetub",
                                   "zar1"),]
gene_pos <- which(rownames(phData) %in% rownames(select))
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,
                                                 labels = rownames(select)))


#接着绘制热图
col_fun = colorRamp2(c(-2, 0, 2), c("cornflowerblue","white","red"))
p1 <- Heatmap(phData,col=col_fun,
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = T,
             show_row_names = FALSE,
             #column_split = levels(cluster_info),
             top_annotation = top_anno,
             right_annotation = row_anno,
             column_title = NULL,
             heatmap_legend_param = list(
               title = "scale",
               title_position = "leftcenter-rot"),column_names_gp = gpar(fontsize = 18),
             row_names_gp = gpar(fontsize = 18))



#phres <- pheatmap(
#  phData, 
#  color = colorRampPalette(c("darkblue", "white", "red3"))(99), #配色
#  scale = "row",
#  cluster_rows = F, #不按行聚类
#  cluster_cols = T, #按列聚类
#  clustering_method = "complete",
#  show_rownames = F, #显示cluster名
#  annotation_row = data.frame(cluster = top50_ovar.new$cluster)), 
#  annotation_colors = list(cluster = cluster_colors)) #注释的配色就用前面设置的cluster的颜色
path <- "/home/tuyx/scRNA-Seq/new_10X/New_0.8/BP_ovary_YDmarked/" ##文件目录
filename <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(filename, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data1 <- lapply(filePath, function(x){
  read_tsv(x,col_names =T)})

dat=purrr::map(1:13,function(i){
  tmp=data1[[i]][which(data1[[i]]$mark==1),];
  })
for (i in 1:13){
  dat[[i]]$cell=gsub(".txt","",names(data1[i]))
}
dat[[13]]
data=rbind(dat[[10]],dat[[5]],dat[[7]],dat[[4]],dat[[13]],
      dat[[11]],dat[[8]],dat[[1]],dat[[3]],dat[[2]],
      dat[[12]],dat[[6]],dat[[9]])
data=data.frame(Description=factor(data$Description,levels = unique(data$Description)),
                qvalue=data$qvalue,
                Count=data$Count,BgRatio=data$BgRatio,
                RichFactor=data$Count/ as.numeric(sub("/\\d+", "",data$BgRatio)),
                CellType=data$cell)
data[,c(1,2,3,5,6)]
p2=ggplot(data[,c(1,2,3,5,6)],aes(CellType,Description)) + geom_point(aes(color=qvalue,size=RichFactor))+
  scale_colour_gradient(low="red",high= "grey")+
  labs(size="RichFactor",x="CellType",y="GO Term",title="Biological processes for ovary")+
  theme_bw()+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(size = 20,angle = 90, vjust = 0.5,hjust = 1))+
  theme(axis.ticks = element_blank())+scale_y_discrete(position = "right")+
  scale_x_discrete(limits=names(col))
ggarrange(p1,p2,align = "v")
