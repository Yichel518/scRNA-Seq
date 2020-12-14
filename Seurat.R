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
tsne_ovary=ovary@reductions$tsne@cell.embeddings
cluster_matrix=as.data.frame(ovary$RNA_snn_res.0.8)
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
ovary <- CreateSeuratObject(counts = ovary_data, project = "ovary", min.cells = 3, 
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
# 选择gene数目小于2500 以及大于200 以及 线粒体数目小于5%的细胞 #可以基于violin图来选择你要卡的阈值
# 过滤完之后还剩余 genes:13714,cells：2638
ovary <- subset(x = ovary, subset = nFeature_RNA > 200 & nFeature_RNA <2000& percent.mito < 20)
#（2）数据标准化
# 使用log转化，度量单位是10000
ovary <- NormalizeData(object = ovary, normalization.method = "LogNormalize", scale.factor = 1e4)

#（3）Identification of highly variable features (feature selection) 识别高度变异基因(特征选择的过程)、
ovary <- FindVariableFeatures(object = ovary,selection.method = 'dispersion', nfeatures = 2000)
# Identify the 10 most highly variable genes 取前十个变化最大的基因
top20 <- head(x = VariableFeatures(object = ovary), 20)
# plot variable features with and without labels 画出2000个高变异基因
plot1 <- VariableFeaturePlot(object = ovary)
LabelPoints(plot = plot1, points = top20, repel = TRUE)

#（4）Scaling 中心化
all.genes <- rownames(x = ovary)
ovary <- ScaleData(object = ovary, features = all.genes)
#- `ScaleData` function to remove unwanted sources of variation from a single-cell dataset, 比如mitochondrial contamination线粒体污染
ovary <- ScaleData(object = ovary, vars.to.regress = 'percent.mito')

# 3.Perform linear dimensional reduction线性降维 + Cluster the cells 对细胞进行聚类分析
##接下来，我们对中心化的数据执行PCA。 默认情况下，只有先前确定的变量要素用作输入，但如果要选择其他子集，则可以使用`features`参数定义。
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
        cols=c("#D2B48C","#DB7093","6A5ACD","#2E8B57","#FFF0F5",
               "#808080","#FFB6C1","#556B2F","#4682B4","#5F9EA0",
               "#DAA520","#BDB76B","#D8BFD8","#F4A460","#A0522D",
               "#C0C0C0","#BA55D3","#4169E1","#008B8B","#FF6347","#800000"))

DimPlot(object = ovary, reduction = 'umap',label = TRUE,pt.size = 1,
        cols=c("#c2ccd0","#808080","#88ada6",
               "#ffb3a7","#4b5cc4","#a1afc9","#a4e2c6","#eedeb0","#ffc773"))


# TSNE
ovary <- RunTSNE(object = ovary, dims = 1:16)
DimPlot(object = ovary, reduction = 'tsne',label = TRUE )
#ovary <- RunICA(object = ovary)
#DimPlot(object = ovary, reduction = 'ica',label = TRUE )
# 您可以在此时保存对象，以便可以轻松地将其重新加载，而无需重新运行上面执行的计算密集型步骤，或者可以轻松地与协作者共享
ovary=readRDS("~/scRNA-Seq/new_10X/ovary_0.8.20201205.rds")
saveRDS(ovary, file = "~/scRNA-Seq/new_10X/ovary_0.8.20201205.rds")


saveRDS(testis, file = "~/scRNA-Seq/Result_testis/testis2.rds")
saveRDS(gonad.integrated, file = "~/scRNA-Seq/gonad_integrated/gonad.integrated2.rds")
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
write.table(ovary_new.markers,"~/scRNA-Seq/new_10X/ovary_new.markers20201205.txt",col.names=T,row.names=T,quote=F,sep="\t")
top50_ovar.new=ovary_new.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) #将ovar.markers传给group_by，group_by按cluster 排序，再将结果传给top_n，top_n按avg_logFC排序，显示每个类中的前两个
write.table(top50_ovar.new,"~/scRNA-Seq/new_10X/top50_ovar.new.txt",col.names=T,row.names=T,quote=F,sep="\t")
top5_ovar.new=ovary_new.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) #将ovar.markers传给group_by，group_by按cluster 排序，再将结果传给top_n，top_n按avg_logFC排序，显示每个类中的前两个
write.table(top5_ovar.new,"~/scRNA-Seq/new_10X/top5_ovar.new.txt",col.names=T,row.names=T,quote=F,sep="\t")

ovary.markers <- FindAllMarkers(object = ovary, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DoHeatmap(object = ovary, features =top5_ovary$gene) + NoLegend()
top5_ovary=ovary_new.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) #将ovar.markers传给group_by，group_by按cluster 排序，再将结果传给top_n，top_n按avg_logFC排序，显示每个类中的前两个
top3_ovary=ovary.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC) #将ovar.markers传给group_by，group_by按cluster 排序，再将结果传给top_n，top_n按avg_logFC排序，显示每个类中的前两个
