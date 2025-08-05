####method-N1####


# library(multtest)
# if(!require(multtest))install.packages("multtest")
# if(!require(Seurat))install.packages("Seurat")
# if(!require(dplyr))install.packages("dplyr")
# if(!require(mindr))install.packages("mindr")
# if(!require(mindr))install.packages("tidyverse")


# #####自动读取cellranger(LINUX)输出的feature barcode matric
# rm(list = ls())
# # 6个样本   Je  Ji  JT  Se  Si  ST
# Je.data <- Read10X(data.dir = "./Je/filtered_feature_bc_matrix") 
# Ji.data <- Read10X(data.dir = "./Ji/filtered_feature_bc_matrix") 
# JT.data <- Read10X(data.dir = "./JT/filtered_feature_bc_matrix") 
# Se.data <- Read10X(data.dir = "./Se/filtered_feature_bc_matrix") 
# Si.data <- Read10X(data.dir = "./Si/filtered_feature_bc_matrix") 
# ST.data <- Read10X(data.dir = "./ST/filtered_feature_bc_matrix") 

# ##### Seurat 4.0 构建Seurat object

# # 6个样本
# Je.obj <- CreateSeuratObject(counts = Je.data, project = "Je")
# Ji.obj <- CreateSeuratObject(counts = Ji.data, project = "Ji")
# JT.obj <- CreateSeuratObject(counts = JT.data, project = "JT")
# Se.obj <- CreateSeuratObject(counts = Se.data, project = "Se")
# Si.obj <- CreateSeuratObject(counts = Si.data, project = "Si")
# ST.obj <- CreateSeuratObject(counts = ST.data, project = "ST")

# # 合并Seurat object
 
# JS.all <- merge(Je.obj, y = c(Ji.obj, JT.obj,Se.obj,Si.obj,ST.obj), add.cell.ids = c("Je", "Ji", "JT",'Se','Si','ST'), project = "J3S3")

####method-N2####
# library(Seurat)
# library(tidyverse)
# library(patchwork)
# # 导入单样本文件
# dir = c(
#   "./Je/filtered_feature_bc_matrix",
#   "./Ji/filtered_feature_bc_matrix",
#   "./JT/filtered_feature_bc_matrix",
#   "./Se/filtered_feature_bc_matrix",
#   "./Si/filtered_feature_bc_matrix",
#   "./ST/filtered_feature_bc_matrix"
# )

# # 读取并创建Seurat Object
# scRNAlist <- list()  # 创建一个空的列表,其中包含三个seurat对象
# names(dir) = c("Je", "Ji", "JT",'Se','Si','ST')
# sample_name = c("Je", "Ji", "JT",'Se','Si','ST')
# for (i in 1:length(dir)){
#   counts <- Read10X(data.dir = dir[i])
#   scRNAlist[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200)
#   scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = sample_name[i])
# }

# # 合并Seurat object
# scRNA <- merge(scRNAlist[[1]],y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]]))
# dim(scRNA)

# # 质量控制
# scRNA <- SCTransform(scRNA, verbose = FALSE)
# scRNA <- RunPCA(scRNA, verbose = FALSE)

# # 去除批次效应
# scRNA <- ScaleData(scRNA, verbose = FALSE)

# # 聚类
# scRNA <- FindNeighbors(scRNA, dims = 1:30)

# scRNA <- FindClusters(scRNA, resolution = 0.5)

# # 降维

# scRNA <- RunUMAP(scRNA, dims = 1:30)

# # 画UMAP
# DimPlot(scRNA, reduction = "umap", label = TRUE) + NoLegend()

# # 画marker genes    

# # 找到marker genes
# default.assay <- DefaultAssay(scRNA)
# markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

####method-N3####
# library(Seurat)
# library(tidyverse)
# library(patchwork)
# # 导入单样本文件
# dir = c(
#   "./Je/filtered_feature_bc_matrix",
#   "./Ji/filtered_feature_bc_matrix",
#   "./JT/filtered_feature_bc_matrix",
#   "./Se/filtered_feature_bc_matrix",
#   "./Si/filtered_feature_bc_matrix",
#   "./ST/filtered_feature_bc_matrix"
# )

# # 读取并创建Seurat Object
# scRNAlist <- list()  # 创建一个空的列表,其中包含三个seurat对象
# names(dir) = c("Je", "Ji", "JT",'Se','Si','ST')
# sample_name = c("Je", "Ji", "JT",'Se','Si','ST')
# for (i in 1:length(dir)){
#   counts <- Read10X(data.dir = dir[i])
#   scRNAlist[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200)
#   scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = sample_name[i])
# }

# # 合并Seurat object
# scRNA <- merge(scRNAlist[[1]],y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]]))
# dim(scRNA)

# #  数据标准化和选择高变基因
# # 每一个样本分别进行数据标准化和提取高变基因
# for (i in 1:length(scRNAlist)){
#   scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
#   scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]],selection.method = 'vst',nfeatures = 2000) # 这里我们只需要2000的基因
# }

# # 以variableFeatures为基础寻找锚点
# scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
# dim(scRNA.anchors)

# # 整合数据
# scRNA.integrated <- IntegrateData(anchorset = scRNA.anchors)
# dim(scRNA.integrated)

# # 缩放
# DefaultAssay(scRNA.integrated) <- "integrated"

# scRNA.integrated <- ScaleData(scRNA.integrated)

# # 聚类

# scRNA.integrated <- RunPCA(scRNA.integrated, verbose = FALSE)

# scRNA.integrated <- RunUMAP(scRNA.integrated, dims = 1:30)

# # 画PCA

# DimPlot(scRNA.integrated, reduction = "pca",group.by = 'orig.ident', label = TRUE) + NoLegend()

# # 肘图，拐点图；# 从图中可以看出我们最好应该选择前n个维度的数据
# ElbowPlot(scRNA.integrated,ndims = 30,reduction = 'pca')

# # 聚类

# scRNA <- FindNeighbors(scRNA.integrated,dims = 1:20)
# scRNA <- FindClusters(scRNA,resolution = 0.1) # 0.5聚类的结果太多
# table(scRNA@meta.data$seurat_clusters)
# metadata <- scRNA@meta.data
# # 单独将数据提取出来
# cell_cluster <-  data.frame(cell_ID=rownames(metadata),cluster_ID=metadata$seurat_clusters)


# # 画UMAP；降维
# scRNA <- RunUMAP(scRNA,dims = 1:20)
# embed_umap <- Embeddings(scRNA,'umap')

# # group_by_cluster
# DimPlot(scRNA,reduction = 'umap')
# # group_by_sample
# DimPlot(scRNA,reduction = 'umap',group.by = 'orig.ident')

# # TSNE降维
# scRNA <- RunTSNE(scRNA,dims = 1:20)
# embed_tsne <- Embeddings(scRNA,'tsne')

# # group_by_cluster
# DimPlot(scRNA,reduction = 'tsne')
# # group_by_sample
# DimPlot(scRNA,reduction = 'tsne',group.by = 'orig.ident')


# # 质控
# # 切换数据集
# DefaultAssay(scRNA) <- 'RNA'
# # 计算线粒体和红细胞的基因比例
# scRNA[['percent.mt']] <- PercentageFeatureSet(scRNA,pattern = 'MT-')

# # 通常是指与血红蛋白(Hemoglobin)相关的基因
# HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
# HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) # 匹配已经拥有的基因,返回一个含有下标的向量
# HB.genes <- rownames(scRNA@assays$RNA)[HB_m]
# scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes)

# ##meta.data添加信息
# proj_name <- data.frame(proj_name=rep("demo2",ncol(scRNA)))
# rownames(proj_name) <- row.names(scRNA@meta.data)
# scRNA <- AddMetaData(scRNA, proj_name)

# # # 画质控指标

# # VlnPlot(scRNA, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'percent.HB'), group.by = 'orig.ident')

# # VlnPlot(scRNA, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'percent.HB'), group.by = 'proj')

# # 可视化
# # 绘制小提琴图
# # 所有样本一个小提琴图用group.by="proj_name"，每个样本一个小提琴图用group.by="orig.ident"
# plot7 <-VlnPlot(scRNA, group.by = "proj_name",  raster=FALSE,
#                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
#                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
# ) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) # 将x轴的标题,文本和刻度线都设置为空,这样x轴就不会显示任何内容
# plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
# pearplot <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none") 
# pearplot

# # 去除细胞特征过高和过低的细胞
# scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

# # 数据归一化
# scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

# # 鉴定高变基因
# # 这一步的目的是鉴定出细胞与细胞之间表达相差很大的基因,用于后续鉴定细胞类型
# # 我们使用默认参数,用vst方法选出2000个高变基因
# scRNA <- FindVariableFeatures(scRNA,selection.method = 'vst',nfeatures = 2000)
# dim(scRNA) # 但是这里跑程序的时候基因的数量不对,还没有找到原因
# # [1] 18037 21402

# # 前十个高变基因
# top10 <- head(VariableFeatures(scRNA),10)
# top10
# #  [1] "PTGDS"    "S100A9"   "S100A8"   "CST3"     "TRBV11-2" "HLA-DQA1" "HLA-DRA"  "C1QA"  "LILRA4"   "LYZ"   

# # 可视化
# plot1 <- VariableFeaturePlot(scRNA)
# plot2 <- LabelPoints(plot = plot1,points = top10,repel = TRUE,xnudge = 0,ynudge = 0)

# # 细胞注释
# Idents(scRNA) <- 'integrated_snn_res.0.5'
# plot3 = DimPlot(scRNA, reduction = "umap", label=T)

# # 鉴定细胞类型
# # 为了后续分析的方便,我们先用singleR来预测每个cluster的细胞类型
# library(celldex)
# library(SingleR)
# cg <- ImmGenData() # 选定一个参考集数据,ImmGenData是一个免疫细胞的数据集
# cellpred <- SingleR(test=testdata,ref=cg, labels=cg$label.main)
# table(cellpred$labels) # 看看都注释到了哪些细胞
# #      B cells      Basophils Endothelial cells       Eosinophils  Epithelial cells 
# #       20337            1               224               595                31 
# #       Fibroblasts               ILC       Macrophages        Mast cells 
# #                3               205                 5                 1 
# # 由得到的结果可以看出,用SingerR注释的结果太拉胯了,虽然手动注释比较麻烦,但是为了数据的准确性还是手动注释吧

# cellType=data.frame(seurat=scRNA@meta.data$seurat_clusters,
#                     predict=cellpred$labels)
# table(cellType[,1:2]) # 访问celltyple的2~3列

# ####细胞注释
# scRNA.sub <- subset(scRNA, idents = c(1,3,4,7), invert = TRUE) # 挑选需要的簇
# scRNA.sub
# new.cluster.ids <- c(
#   "1" = "CD4 T",
#   "2" = "Monocyte",
#   "3" = "CD4 T",
#   "4" = "NK",
#   "5" = "Monocyte",
#   "6" = "CD8 T",
#   "7" = "B",
#   "8" = "CD8 T",
#   "9" ="Monocyte",
#   "11"="NK",
#   "13" = "Monocyte",
#   "14" = "Monocyte",
#   "16" = "CD4 T",
#   "17" = "DC",
#   "18" = "B",
#   "19"="NK",
#   "20"="Monocyte",
#   "21"="DC",
#   "22" = "B",
#   "25"="NK"
# )
# scRNA.sub <- RenameIdents(scRNA.sub, new.cluster.ids)
# options(repr.plot.height = 6, repr.plot.width = 8)
# scRNA.sub$cell_type <- Idents(scRNA.sub)
# Idents(scRNA.sub) <- "cell_type"
# DimPlot(scRNA.sub, reduction = "umap", label = TRUE,raster=FALSE)

# # 保存图片
# ggsave(filename = "../output/images/scRNA3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)
# # 保存数据
# saveRDS(scRNA, file = "../output/scRNA3k_final.rds")

####method-N4####
# library(Seurat)
# library(tidyverse)
# library(patchwork)
# # 导入单样本文件
# dir = c(
#   "./Je/filtered_feature_bc_matrix",
#   "./Ji/filtered_feature_bc_matrix",
#   "./JT/filtered_feature_bc_matrix",
#   "./Se/filtered_feature_bc_matrix",
#   "./Si/filtered_feature_bc_matrix",
#   "./ST/filtered_feature_bc_matrix"
# )

# # 读取并创建Seurat Object
# scRNAlist <- list()  # 创建一个空的列表,其中包含三个seurat对象
# names(dir) = c("Je", "Ji", "JT",'Se','Si','ST')
# sample_name = c("Je", "Ji", "JT",'Se','Si','ST')
# for (i in 1:length(dir)){
#   counts <- Read10X(data.dir = dir[i])
#   scRNAlist[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200)
#   scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = sample_name[i])
# }

# # 合并Seurat object
# scRNA <- merge(scRNAlist[[1]],y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]]))
# dim(scRNA)

# # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

# # Visualize QC metrics

# VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# # FeatureScatter

# plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# # Combine the plots and print them all together
# CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none")

# # 去除细胞特征过高和过低的细胞

# scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

# # 数据归一化

# scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

# # 鉴定高变基因

# # 这一步的目的是鉴定出细胞与细胞之间表达相差很大的基因,用于后续鉴定细胞类型

# # 我们使用默认参数,用vst方法选出2000个高变基因

# scRNA <- FindVariableFeatures(scRNA, selection.method = 'vst', nfeatures = 2000)

# # 前十个高变基因

# top10 <- head(VariableFeatures(scRNA),10)
# top10

# # 可视化

# plot1 <- VariableFeaturePlot(scRNA)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

# # 细胞注释

# Idents(scRNA) <- 'SNN_res.0.5'

# plot3 = DimPlot(scRNA, reduction = "umap", label=T)

# # 鉴定细胞类型

# # 为了后续分析的方便, 我们先用singleR来预测每个cluster的细胞类型

# library(celldex)
# library(SingleR)

# cg <- ImmGenData() # 选定一个参考集数据, ImmGenData是一个免疫细胞的数据集

# cellpred <- SingleR(test=scRNA,ref=cg, labels=cg$label.main)

# table(cellpred$labels) # 看看都注释到了哪些细胞

# #      B cells      Basophils Endothelial cells       Eosinophils  Epithelial cells

# #       20337            1               224               595                31

# #       Fibroblasts               ILC       Macrophages        Mast cells

# #                3               205                 5                 1

# # 由得到的结果可以看出,用SingerR注释的结果太拉胯了,虽然手动注释比较麻烦,但是为了数据的准确性还是手动注释吧

# cellType=data.frame(seurat=Idents(scRNA),

#                      predict=cellpred$labels)

# table(cellType[,1:2]) # 访问celltyple的2~3列

# ####细胞注释

# scRNA.sub <- subset(scRNA, idents = c(1,3,4,7), invert = TRUE) # 挑选需要的簇

# scRNA.sub

# new.cluster.ids <- c(
#     "1" = "CD4 T",
#     "2" = "Monocyte",
#     "3" = "CD4 T",
#     "4" = "NK",
#     "5" = "Monocyte",
#     "6" = "CD8 T",
#     "7" = "B",
#     "8" = "CD8 T",
#     "9" ="Monocyte",
#     "11"="NK",
#     "13" = "Monocyte",
#     "14" = "Monocyte",
#     "16" = "CD4 T",
#     "17" = "DC",
#     "18" = "B",
#     "19"="NK",
#     "20"="Monocyte",
#     "21"="DC",
#     "22" = "B",
#     "25"="NK"
# )

# scRNA.sub <- RenameIdents(scRNA.sub, new.cluster.ids)

# options(repr.plot.height = 6, repr.plot.width = 8)

# scRNA.sub$cell_type <- Idents(scRNA.sub)

# Idents(scRNA.sub) <- "cell_type"

# DimPlot(scRNA.sub, reduction = "umap", label = TRUE,raster=FALSE)

# # 保存图片

# ggsave(filename = "../output/images/scRNA3k_umap_subcluster.jpg", height = 7, width = 12, plot = plot, quality = 50)

# # 保存数据

# saveRDS(scRNA.sub, file = "../output/scRNA3k_subcluster_final.rds")
# 小鼠线粒体基因：^mt-
####method-N5####
library(Seurat)
library(tidyverse)
library(patchwork)
# 导入单样本文件
dir = c(
  "./Je/filtered_feature_bc_matrix",
  "./Ji/filtered_feature_bc_matrix",
  "./JT/filtered_feature_bc_matrix",
  "./Se/filtered_feature_bc_matrix",
  "./Si/filtered_feature_bc_matrix",
  "./ST/filtered_feature_bc_matrix"
)

# 读取并创建Seurat Object
# S手术组，J假手术组
scRNAlist <- list()  # 创建一个空的列表,其中包含三个seurat对象
names(dir) = c("Je", "Ji", "JT",'Se','Si','ST')
sample_name = c("Je", "Ji", "JT",'Se','Si','ST')
for (i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = sample_name[i])
}

# 合并Seurat object
scRNA <- merge(scRNAlist[[1]],y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]]))
dim(scRNA)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
# Visualize QC metrics as a violin plot
p_QC_n1 <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "p_QC_n1.pdf", height = 7, width = 12, plot = p_QC_n1)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p_QC_FeatureScatter_n2 <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none")
ggsave(filename = "p_QC_FeatureScatter_n2.pdf", height = 7, width = 12, plot = p_QC_FeatureScatter_n2)

scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)

# Normalization
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scRNA), 10)
plot3 <- VariableFeaturePlot(scRNA)

# Label features of interest
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
p_ideFeatures <- CombinePlots(plots = list(plot3, plot4), nrow=1, legend="none")
ggsave(filename = "p_ideFeatures.pdf", height = 7, width = 12, plot = p_ideFeatures)

# Scale data
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA)

# Run PCA
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))

# Examine PCA results
print(scRNA[["pca"]], dims = 1:5, nfeatures = 5)
p_PCA_n1 <- VizDimLoadings(scRNA, dims = 1:2, reduction = "pca")
ggsave(filename = "p_PCA_n1.pdf", height = 7, width = 12, plot = p_PCA_n1)
p_PCA_n2 <- DimPlot(scRNA, reduction = "pca")
ggsave(filename = "p_PCA_n2.pdf", height = 7, width = 12, plot = p_PCA_n2)

pdf(file="p_Heatmap1.pdf",,height = 7, width = 12)
DimHeatmap(scRNA, dims = 1, cells = 500, balanced = TRUE)
dev.off()

pdf(file="p_Heatmap15.pdf",,height = 7, width = 12)
DimHeatmap(scRNA, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# Determine the optimal number of principal components
scRNA <- JackStraw(scRNA, num.replicate = 100)

p_Elbow <- ElbowPlot(scRNA)
ggsave(filename = "p_Elbow.pdf", height = 7, width = 12, plot = p_Elbow)

# Cluster the cells
scRNA <- FindNeighbors(scRNA, dims = 1:15)
scRNA <- FindClusters(scRNA, resolution = 0.5)

# Run UMAP/TSNE
scRNA <- RunUMAP(scRNA, dims = 1:15)
# Clustering results can be visualized using either UMAP or tSNE
p_UMAP <- DimPlot(scRNA, reduction = "umap", label = TRUE, raster=FALSE)
ggsave(filename = "p_UMAP.pdf", height = 7, width = 12, plot = p_UMAP)

p_UMAP_orig.ident <- DimPlot(scRNA, reduction = "umap", group.by='orig.ident')
ggsave(filename = "p_UMAP_orig.ident.pdf", height = 7, width = 12, plot = p_UMAP_orig.ident)

# Save the Seurat object
# saveRDS(scRNA, file = "JS_seurat_tutorial.rds")

# # find all markers of cluster 2
# cluster2.markers <- FindMarkers(scRNA, ident.1 = 2)
# head(cluster2.markers, n = 5)

# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(scRNA, ident.1 = 5, ident.2 = c(0, 3))
# head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE)
# saveRDS(scRNA, file = "JS_scRNA.markers.rds")

scRNA.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

# # cluster0.markers <- FindMarkers(scRNA, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
# p_markers_vln1 <- VlnPlot(scRNA, features = c("MS4A1", "CD79A"))
# ggsave(filename = "p_markers_vln1.pdf", height = 7, width = 12, plot = p_markers_vln1)

# # you can plot raw counts as well
# VlnPlot(scRNA, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

# FeaturePlot(scRNA, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
#     "CD8A"))

scRNA.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
pdf(file="p_clustertop10gene.pdf",,height = 20, width = 30)
DoHeatmap(scRNA, features = top10$gene) + NoLegend()
dev.off()
# saveRDS(scRNA, file = "JS_seurat_cluster_filt.rds")
# 细胞注释
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
DimPlot(scRNA, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(scRNA, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "../output/images/scRNA3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

# saveRDS(scRNA, file = "JS_seurat_final.rds")



# 展示细胞簇的基因表达情况

genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','TNFRSF18','SLC4A10','IL7R',
                   'CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',
                   'CD68', 'CD163', 'CD14', 'JCHAIN',
                   'TPSAB1' , 'TPSB2',
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'C1QA',  'C1QB', 
                   'S100A9', 'S100A8','CD86', 'MMP19',
                   'LAMP3', 'IDO1','IDO2',
                   'CD1E','CD1C',
                   'KRT86','GNLY',
                   'FGF7','MME','ACTA2','GFPT2','CNN1','CNN2',
                   'DCN', 'LUM',  'GSN' , 
                   'FAP','FN1','THY1','COL1A1','COL3A1', 
                   'PECAM1', 'VWF',
                   'EPCAM' , 'KRT19', 'PROM1', 'CD24','MKI67',
                   'ALDH1A1', 'ALDH1A3','ITGA4','ITGA6')
pdf(file="p_clustergenedotplot.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(genes_to_check),cluster.idents = T) + coord_flip()
dev.off()


# 根据基因表达图手动注释细胞
B_P_cell=c(1,25,36)
Mast_cell=c(39)
T_cell=c(0,2)
Epithelial=c(14,10,11,23,31,46,21,28,30,35,32,7,18,16,26,17,19,43,45,22,12,33,42,20,44,8,29,6,9,4,15,40,13,37)
Endothelial =c(24)
Myeloid=c(3)
Fibro_cell=c(5,47)
MIX=c(38,27,34,41)
# table(c(B_P_cell,Mast_cell,T_cell,Epithelial,Endothelial,Myeloid,Fibro_cell,MIX))
# print(setdiff(0:47,c(B_P_cell,Mast_cell,T_cell,Epithelial,Endothelial,Myeloid,Fibro_cell,MIX)))
current.cluster.ids <- c(B_P_cell,Mast_cell,T_cell,Epithelial,Endothelial,Myeloid,Fibro_cell,MIX)
new.cluster.ids <- c(rep("B cell",length(B_P_cell)),
                     rep("Mast cell",length(Mast_cell)),
                     rep("T cell",length(T_cell)),
                     rep("Epithelial cell",length(Epithelial)),
                     rep("Endothelial cell",length(Endothelial)),
                     rep("Myeloid cell",length(Myeloid)),
                     rep("Fibrocyte",length(Fibro_cell)),
                     rep("MIX",length(MIX))
                     )
scRNA@meta.data$pre_cellType <- plyr::mapvalues(x = scRNA$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
scRNA <- subset(scRNA, pre_cellType %in% c("T cell","B cell","Myeloid cell","Epithelial cell","Fibrocyte",
                                                               "Endothelial cell","Mast cell"))


# 合并数据自动注释
#合并数据，并进行预处理
BRCA_SingleCell.T <- merge(G1, G2)
rm(G1)
rm(G2)
BRCA_SingleCell.T <- NormalizeData(BRCA_SingleCell.T, normalization.method = "LogNormalize")
BRCA_SingleCell.T <- FindVariableFeatures(BRCA_SingleCell.T, selection.method = "vst", nfeatures = 3000)
BRCA_SingleCell.T <- ScaleData(BRCA_SingleCell.T, vars.to.regress = c('nCount_RNA'))
BRCA_SingleCell.T <- RunPCA(BRCA_SingleCell.T)
BRCA_SingleCell.T <- RunHarmony(object = BRCA_SingleCell.T,group.by.vars = c('Dataset'))
BRCA_SingleCell.T <- RunUMAP(BRCA_SingleCell.T,reduction = "harmony",dims = 1:30,seed.use = 12345)
BRCA_SingleCell.T <- FindNeighbors(BRCA_SingleCell.T,reduction = 'harmony', dims = 1:30, verbose = FALSE)
BRCA_SingleCell.T <- FindClusters(BRCA_SingleCell.T,resolution = 0.5, verbose = FALSE,random.seed=20220727)
#对合并数据重新进行细胞自动注释
cellAnnotation <- function(obj,markerList,assay='SCT',slot='data'){
    markergene <- unique(do.call(c,markerList))
    Idents(obj) <- 'seurat_clusters'
    cluster.averages <- AverageExpression(obj,assays=assay, slot = slot,features=markergene, return.seurat = TRUE)
    if(assay=='SCT'){
        scale.data <- cluster.averages@assays$SCT@data
    }else{
        scale.data <- cluster.averages@assays$RNA@data
    }
    print(scale.data)
    cell_score <- sapply(names(markergeneList),function(x){
        tmp <- scale.data[rownames(scale.data)%in%markergeneList[[x]],]
        if(is.matrix(tmp)){
            if(nrow(tmp)>=2){
                res <- apply(tmp,2,max)
                return(res)
            }else{
                return(rep(-2,ncol(tmp)))
            }
        }else{
            return(tmp)
        }
    })
    print(cell_score)
    celltypeMap <- apply(cell_score,1,function(x){
        colnames(cell_score)[which(x==max(x))]
    },simplify = T)
    obj@meta.data$cellType_auto <- plyr::mapvalues(x = obj@active.ident, from = names(celltypeMap), to = celltypeMap)
    return(obj)
}
lymphocyte <- c('CD3D','CD3E','CD79A','MS4A1','MZB1')
myeloid <- c('CD68','CD14','TPSAB1' , 'TPSB2','CD1E','CD1C','LAMP3', 'IDO1')
EOC <- c('EPCAM','KRT19','CD24')
fibo_gene <- c('DCN','FAP','COL1A2')
endo_gene <- c('PECAM1','VWF')
markergeneList <- list(lymphocyte=lymphocyte,myeloid=myeloid,Epi=EOC,fibo=fibo_gene,endo=endo_gene)
BRCA_SingleCell.T <- cellAnnotation(obj=BRCA_SingleCell.T,assay = 'RNA',markerList=markergeneList)                                                               



####method-N6####
library(Seurat)
library(tidyverse)
library(patchwork)
# 导入单样本文件
dir = c(
  "./Je/filtered_feature_bc_matrix",
  "./Ji/filtered_feature_bc_matrix",
  "./JT/filtered_feature_bc_matrix",
  "./Se/filtered_feature_bc_matrix",
  "./Si/filtered_feature_bc_matrix",
  "./ST/filtered_feature_bc_matrix"
)

# 读取并创建Seurat Object
# S手术组，J假手术组
scRNAlist <- list()  # 创建一个空的列表,其中包含三个seurat对象
names(dir) = c("Je", "Ji", "JT",'Se','Si','ST')
sample_name = c("Je", "Ji", "JT",'Se','Si','ST')
for (i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = sample_name[i])
}

# 合并Seurat object
scRNA <- merge(scRNAlist[[1]],y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]]))
dim(scRNA)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
# Visualize QC metrics as a violin plot
p_QC_n1 <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "p_QC_n1.pdf", height = 7, width = 12, plot = p_QC_n1)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p_QC_FeatureScatter_n2 <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none")
ggsave(filename = "p_QC_FeatureScatter_n2.pdf", height = 7, width = 12, plot = p_QC_FeatureScatter_n2)

scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)

# Normalization
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

scRNA <- FindVariableFeatures(scRNA,nfeatures = 3000)
scRNA <- ScaleData(scRNA,vars.to.regress = 'percent.mt')
scRNA <- RunPCA(scRNA)
scRNA <- FindNeighbors(scRNA, reduction = "pca", dims = 1:10)
scRNA <- FindClusters(scRNA, resolution = 0.3)
scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:10) 
pdf(file="p_DimPlotUMAPBH1.pdf",,height = 7, width = 12)
DimPlot(scRNA, reduction = "umap", label = T)
dev.off()
pdf(file="p_DimPlotUMAPBH2.pdf",,height = 7, width = 12)
DimPlot(scRNA,group.by = 'orig.ident')
dev.off()
head(scRNA@meta.data)
pdf(file="p_FeaturePlotUMAPBH.pdf",,height = 7, width = 12)
FeaturePlot(scRNA,features = c('nCount_Spatial', 'nFeature_Spatial', 'percent.mt'))
dev.off()

scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident")
scRNA  <- FindNeighbors(scRNA , dims = 1:30 , reduction = "harmony")
scRNA  <- FindClusters(scRNA, save.snn=T , resolution = 0.3)
scRNA  <- RunUMAP(scRNA , dims=1:30,reduction='harmony')

pdf(file="p_DimPlotUMAPAH1.pdf",height = 7, width = 12)
DimPlot(scRNA , reduction = "umap",group.by = 'seurat_clusters',label = T)
dev.off()
pdf(file="p_DimPlotUMAPAH2.pdf",,height = 7, width = 12)
DimPlot(scRNA,group.by = 'orig.ident')
dev.off()

# # Identify highly variable features
# scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
# top10 <- head(VariableFeatures(scRNA), 10)
# plot3 <- VariableFeaturePlot(scRNA)

# # Label features of interest
# plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
# p_ideFeatures <- CombinePlots(plots = list(plot3, plot4), nrow=1, legend="none")
# ggsave(filename = "p_ideFeatures.pdf", height = 7, width = 12, plot = p_ideFeatures)

# # Scale data
# all.genes <- rownames(scRNA)
# scRNA <- ScaleData(scRNA)

# # Run PCA
# scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))

# # Examine PCA results
# print(scRNA[["pca"]], dims = 1:5, nfeatures = 5)
# p_PCA_n1 <- VizDimLoadings(scRNA, dims = 1:2, reduction = "pca")
# ggsave(filename = "p_PCA_n1.pdf", height = 7, width = 12, plot = p_PCA_n1)
# p_PCA_n2 <- DimPlot(scRNA, reduction = "pca")
# ggsave(filename = "p_PCA_n2.pdf", height = 7, width = 12, plot = p_PCA_n2)

# pdf(file="p_Heatmap1.pdf",,height = 7, width = 12)
# DimHeatmap(scRNA, dims = 1, cells = 500, balanced = TRUE)
# dev.off()

# pdf(file="p_Heatmap10.pdf",,height = 7, width = 12)
# DimHeatmap(scRNA, dims = 1:10, cells = 500, balanced = TRUE)
# dev.off()

# # Determine the optimal number of principal components
# scRNA <- JackStraw(scRNA, num.replicate = 100)

# p_Elbow <- ElbowPlot(scRNA)
# ggsave(filename = "p_Elbow.pdf", height = 7, width = 12, plot = p_Elbow)

# # Cluster the cells
# scRNA <- FindNeighbors(scRNA, dims = 1:10)
# scRNA <- FindClusters(scRNA, resolution = 0.5)

# # Run UMAP/TSNE
# scRNA <- RunUMAP(scRNA, dims = 1:10)
# # Clustering results can be visualized using either UMAP or tSNE
# p_UMAP <- DimPlot(scRNA, reduction = "umap", label = TRUE, raster=FALSE)
# ggsave(filename = "p_UMAP.pdf", height = 7, width = 12, plot = p_UMAP)

# p_UMAP_orig.ident <- DimPlot(scRNA, reduction = "umap", group.by='orig.ident')
# ggsave(filename = "p_UMAP_orig.ident.pdf", height = 7, width = 12, plot = p_UMAP_orig.ident)

# Save the Seurat object
# saveRDS(scRNA, file = "JS_seurat_tutorial_AH.rds")

# # find all markers of cluster 2
# cluster2.markers <- FindMarkers(scRNA, ident.1 = 2)
# head(cluster2.markers, n = 5)

# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(scRNA, ident.1 = 5, ident.2 = c(0, 3))
# head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE)
scRNA.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

# # cluster0.markers <- FindMarkers(scRNA, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
# p_markers_vln1 <- VlnPlot(scRNA, features = c("MS4A1", "CD79A"))
# ggsave(filename = "p_markers_vln1.pdf", height = 7, width = 12, plot = p_markers_vln1)

# # you can plot raw counts as well
# VlnPlot(scRNA, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

# FeaturePlot(scRNA, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
#     "CD8A"))

scRNA.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
pdf(file="p_clustertop10geneAH.pdf",,height = 20, width = 30)
DoHeatmap(scRNA, features = top10$gene) + NoLegend()
dev.off()
# saveRDS(scRNA, file = "JS_seurat_cluster_filt_AH.rds")

features<-c('Adipoq','Pdgfra','Msin','Jam2','Prox1','Steap4','Myocd','Mafb','Cybb','FIt3','Cpa3','Csf3r','Ms4a1','KIrd1','II7r','Dcdc2a','Erbb4')
pdf(file="p_FeaturePlot.pdf",,height = 20, width = 30)
FeaturePlot(object = scRNA, features = features)
dev.off()


####method-N7####
library(Seurat)
library(tidyverse)
library(patchwork)
# 导入单样本文件
dir = c(
  "./Je/filtered_feature_bc_matrix",
  "./Ji/filtered_feature_bc_matrix",
  "./JT/filtered_feature_bc_matrix",
  "./Se/filtered_feature_bc_matrix",
  "./Si/filtered_feature_bc_matrix",
  "./ST/filtered_feature_bc_matrix"
)

# 读取并创建Seurat Object
scRNAlist <- list()  # 创建一个空的列表,其中包含三个seurat对象
names(dir) = c("Je", "Ji", "JT",'Se','Si','ST')
sample_name = c("Je", "Ji", "JT",'Se','Si','ST')
for (i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = sample_name[i])
}

# 合并Seurat object
# scRNA <- merge(scRNAlist[[1]],y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]]))
# dim(scRNA)

#  数据标准化和选择高变基因
# 每一个样本分别进行数据标准化和提取高变基因
features <- SelectIntegrationFeatures(object.list = scRNAlist)
for (i in 1:length(scRNAlist)){
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-")
    scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
    scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
    scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]],selection.method = 'vst',nfeatures = 3000) # 这里我们只需要2000的基因
    scRNAlist[[i]] <- ScaleData(scRNAlist[[i]], features = features, verbose = FALSE)
    scRNAlist[[i]] <- RunPCA(scRNAlist[[i]], features = features, verbose = FALSE)
}

scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, dims = 1:30)
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)

DefaultAssay(scRNA) <- "integrated"

# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, verbose = FALSE)
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
pdf(file="p_ElbowPlot_Anchors.pdf",,height = 7, width = 12)
ElbowPlot(scRNA,ndims = 30)
dev.off()
# t-SNE and Clustering
# scRNA <- RunTSNE(scRNA, reduction = "pca", dims = 1:30)
scRNA  <- RunUMAP(scRNA , dims=1:30,reduction='pca')
scRNA <- FindNeighbors(scRNA, reduction = "pca", dims = 1:30)
scRNA <- FindClusters(scRNA,resolution = 0.8)
saveRDS(scRNA,"scRNA_rmMT_update.seurat.RDS")

scRNA <- readRDS("scRNA_rmMT_update.seurat.RDS")
# 查看各样本信息
table(scRNA@meta.data$orig.ident)
p1=DimPlot(scRNA,group.by = "orig.ident")
# 查看聚类信息
p2=DimPlot(scRNA,group.by = "seurat_clusters",label=T)
combine<-CombinePlots(list(p1,p2))
ggsave(filename = "p_DimPlot_sample_seurat_clusters.pdf", plot = combine, height = 10, width = 14)
#####按样本拆分scRNA#####
scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))
p1=DimPlot(scRNASeJe,group.by = "orig.ident")
# 查看聚类信息
p2=DimPlot(scRNASeJe,group.by = "seurat_clusters",label=T)
combine<-CombinePlots(list(p1,p2))
ggsave(filename = "p_DimPlot_sample_seurat_clusters-scRNASeJe.pdf", plot = combine, height = 10, width = 14)

scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
p1=DimPlot(scRNASiJi,group.by = "orig.ident")
# 查看聚类信息
p2=DimPlot(scRNASiJi,group.by = "seurat_clusters",label=T)
combine<-CombinePlots(list(p1,p2))
ggsave(filename = "p_DimPlot_sample_seurat_clusters-scRNASiJi.pdf", plot = combine, height = 10, width = 14)

scRNASTJT=subset(x = scRNA, subset = (orig.ident == "ST"|orig.ident == "JT"))
p1=DimPlot(scRNASTJT,group.by = "orig.ident")
# 查看聚类信息
p2=DimPlot(scRNASTJT,group.by = "seurat_clusters",label=T)
combine<-CombinePlots(list(p1,p2))
ggsave(filename = "p_DimPlot_sample_seurat_clusters-scRNASTJT.pdf", plot = combine, height = 10, width = 14)

# 根据基因表达图手动注释细胞
Adipocyte_Adipose_Tissue=c(0,1,2,7,8,27)
# Astrocyte=c(27)
T_cell=c(11)
B_cell=c(10,28)
Endothelial_Cell =c(4,5,24,26,31)
Fibroblast=c(6,13,14,25)
Macrophage=c(3,18,19,34,35)
# Microglial_Cell=c(35)
Dendritic_Cell=c(9,15,16,33,36)
Muscle_Cell=c(12,20,29)
Vascular_Stem_Cell=c(23)
Neuron_Cell=c(30,32)
Smooth_Muscle_Cell=c(17,21)
NK_Cell=c(22)

# table(c(B_P_cell,Mast_cell,T_cell,Epithelial,Endothelial,Myeloid,Fibro_cell,MIX))
# print(setdiff(0:47,c(B_P_cell,Mast_cell,T_cell,Epithelial,Endothelial,Myeloid,Fibro_cell,MIX)))

current.cluster.ids <- c(Adipocyte_Adipose_Tissue,T_cell,B_cell,Endothelial_Cell,Fibroblast,Macrophage,Dendritic_Cell,Muscle_Cell,Vascular_Stem_Cell,Neuron_Cell,Smooth_Muscle_Cell,NK_Cell)
new.cluster.ids <- c(rep("Adipocyte Adipose Tissue",length(Adipocyte_Adipose_Tissue)),
                     rep("T cell",length(T_cell)),
                     rep("B cell",length(B_cell)),
                     rep("Endothelial cell",length(Endothelial_Cell)),
                     rep("Fibroblast",length(Fibroblast)),
                     rep("Macrophagel",length(Macrophage)),
                     rep("Dendritic Cell",length(Dendritic_Cell)),
                     rep("Muscle Cell",length(Muscle_Cell)),
                     rep("Vascular Stem Cell",length(Vascular_Stem_Cell)),
                     rep("Neuron Cell",length(Neuron_Cell)),
                     rep("Smooth Muscle Cell",length(Smooth_Muscle_Cell)),
                     rep("NK Cell",length(NK_Cell))
                     )
scRNA@meta.data$cellType <- plyr::mapvalues(x = scRNA$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
# scRNA <- subset(scRNA, pre_cellType %in% c("T cell","B cell","Myeloid cell","Epithelial cell","Fibrocyte","Endothelial cell","Mast cell",'MIX'))

pdf(file="p_cellType_anno.pdf",height = 15, width = 15)
DimPlot(scRNA,group.by = "cellType",label=T)
dev.off()

scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE)
# scRNA.markers %>%
#     group_by(cluster) %>%
# #     dplyr::filter(avg_log2FC > 1)
# saveRDS(scRNA.markers, file = "JS_scRNA.markers.rds")
# scRNA.markers %>%
#     group_by(cluster) %>%
#     dplyr::filter(avg_log2FC > 1) %>%
#     slice_head(n = 50) %>%
#     ungroup() -> top50

# pdf(file="p_clustertop50geneAnchors.pdf",,height = 20, width = 30)
# DoHeatmap(scRNA, features = top50$gene) + NoLegend()
# dev.off()
Adipocyte_Adipose_Tissue=c(0,1,2,7,8,27)
# Astrocyte=c(27)
T_cell=c(11)
B_cell=c(10,28)
Endothelial_Cell =c(4,5,24,26,31)
Fibroblast=c(6,13,14,25)
Macrophage=c(3,18,19,34,35)
# Microglial_Cell=c(35)
Dendritic_Cell=c(9,15,16,33,36)
Muscle_Cell=c(12,20,29)
Vascular_Stem_Cell=c(23)
Neuron_Cell=c(30,32)
Smooth_Muscle_Cell=c(17,21)
NK_Cell=c(22)
clusterAdipocyte.markers <- FindMarkers(scRNA, ident.1 = c(0,1,2,7,8,27),only.pos = TRUE)
clusterT_cell.markers <- FindMarkers(scRNA, ident.1 = c(11),only.pos = TRUE)
clusterB_cell.markers <- FindMarkers(scRNA, ident.1 = c(10,28),only.pos = TRUE)
clusterEndothelia.markers <- FindMarkers(scRNA, ident.1 = c(4,5,24,26,31),only.pos = TRUE)
clusterFibroblast.markers <- FindMarkers(scRNA, ident.1 = c(6,13,14,25),only.pos = TRUE)
clusterMacrophage.markers <- FindMarkers(scRNA, ident.1 = c(3,18,19,34,35),only.pos = TRUE)
clusterDendritic.markers <- FindMarkers(scRNA, ident.1 = c(9,15,16,33,36),only.pos = TRUE)
clusterMuscle.markers <- FindMarkers(scRNA, ident.1 = c(12,20,29),only.pos = TRUE)
clusterVascular_Stem.markers <- FindMarkers(scRNA, ident.1 = c(23),only.pos = TRUE)
clusterNeuron.markers <- FindMarkers(scRNA, ident.1 = c(30,32),only.pos = TRUE)
clusterSmooth_Muscle.markers <- FindMarkers(scRNA, ident.1 = c(17,21),only.pos = TRUE)
clusterNK_Cell.markers <- FindMarkers(scRNA, ident.1 = c(22),only.pos = TRUE)
write.csv(clusterAdipocyte.markers,'clusterAdipocyte.markers.csv')
write.csv(clusterT_cell.markers,'clusterT_cell.markers.csv')
write.csv(clusterB_cell.markers,'clusterB_cell.markers.csv')
write.csv(clusterEndothelia.markers,'clusterEndothelia.markers.csv')
write.csv(clusterFibroblast.markers,'clusterFibroblast.markers.csv')
write.csv(clusterMacrophage.markers,'clusterMacrophage.markers.csv')
write.csv(clusterDendritic.markers,'clusterDendritic.markers.csv')
write.csv(clusterMuscle.markers,'clusterMuscle.markers.csv')
write.csv(clusterVascular_Stem.markers,'clusterVascular_Stem.markers.csv')
write.csv(clusterNeuron.markers,'clusterNeuron.markers.csv')
write.csv(clusterSmooth_Muscle.markers,'clusterSmooth_Muscle.markers.csv')
write.csv(clusterNK_Cell.markers,'clusterNK_Cell.markers.csv')

# 20241227
genes_to_check = c("Nrg4","Nnat","Ghr","Tenm4","Adrb3",
                  "Bank1","Ralgps2","Igkc","Ighm","Bcl11a",
                  "Plxnc1","Gramd1b","Dscam","Samsn1","Tmtc2",
                  "Cyyr1","Cdh13","Ptprb","Etl4","Dach1",
                  "Dcn","Tmeff2","Lama2","Egfr","Col3a1",
                  "F13a1","Mrc1","Rbpj","Mctp1","Dab2",
                  "Ttn","Ctnna3","Trdn","Mylk4","Neb",
                  "Gria4","Ntng1","Klf5","Pde1c","Map2",
                  "Gm2682","Klra9","Klrk1","Klre1","Ptprc",
                  "Cacna1c","Prkg1","Notch3","Myh11","Kcnq5",
                  "Skap1","Gm2682","Ms4a4b","Themis","Grap2",
                  "Top2a","Kif4","Mki67","Prc1","Kif23")


scRNA@meta.data$cellType=factor(scRNA@meta.data$cellType,levels =c('Adipocyte Adipose Tissue','B cell','Dendritic Cell','Endothelial cell','Fibroblast','Macrophagel','Muscle Cell','Neuron Cell','NK Cell','Smooth Muscle Cell','T cell','Vascular Stem Cell'))
pdf(file="p_celltypegenedotplot.pdf",height = 15, width = 15)
DotPlot(scRNA,group.by = 'cellType', features = unique(genes_to_check)) + coord_flip()
dev.off()
pdf(file="p_clustersgenedotplot.pdf",height = 15, width = 15)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(genes_to_check)) + coord_flip()
dev.off()

scRNA <- readRDS("scRNA_rmMT_update.seurat.RDS")

# 20250103
genes_to_check = c("Adipoq","Pnpla2","Plin1","Cidec","Apoc1","Fabp4",
                  "Pxk","Ms4a1","Cd19","Cd74","Cd79a","Ighd",
                  "Cd83","Cd86","Ly75",
                  "Itgax","Zbtb46","Cd86","Cd83","Cd1a",
                  "Cd93","Vwf","Emcn","Egfl7","Flt1","Id3",
                  "Vim","Pdgfrb","Lum","Col6a2","Vtn","Mfap5",
                  "Cd68","Fcgr1","Lyz2","Sepp1","Naaa","Ccr2","Cd74",
                  "Tnnt3","Ttn",
                  "Map2",
                  "Nkg7","Klrf1","Klrd1","Gnly","Ncr1",
                  "Acta2","Myl9","Rgs5","Mylk","Nebl","Myh11",
                  "Trbc2","Cd3d","Cd3g","Cd3e","Il7r","Ltb",
                  "Cd34","Cd31","Cd133","Vegfr2","Vwf")

scRNA@meta.data$seurat_clusters=factor(scRNA@meta.data$seurat_clusters,levels =c("0","1","2","7","8","27","10","28","4","5","24","26","31","6","13","14","25","3","18","19","33","34","12","29","22","17","21","9","15","16","20","11","30","32","35","36","23"))
pdf(file="p_clustersgenedotplot2.pdf",height = 15, width = 15)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(genes_to_check)) + coord_flip()
dev.off()
pdf(file="p_clustersgenedotplot2.pdf",height = 15, width = 15)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(genes_to_check)) + coord_flip()
dev.off()

# 根据基因表达图手动注释细胞
Adipocyte_Adipose_Tissue=c(0,1,2,7,8,27)
# Astrocyte=c(27)
B_cell=c(10,28)
Endothelial_Cell =c(4,5,24,26,31)
Fibroblast=c(6,13,14,25)
Macrophage=c(33)
Monocyte_Cell =c(34)
MP=c(3,18,19)
Muscle_Cell=c(12,29)
myoFB=c(9,15,16,20)
NK_Cell=c(22)
Smooth_Muscle_Cell=c(17,21)
T_cell=c(11)
unknown=c(30,32,35,36,23)

# table(c(B_P_cell,Mast_cell,T_cell,Epithelial,Endothelial,Myeloid,Fibro_cell,MIX))
# print(setdiff(0:47,c(B_P_cell,Mast_cell,T_cell,Epithelial,Endothelial,Myeloid,Fibro_cell,MIX)))

current.cluster.ids <- c(Adipocyte_Adipose_Tissue,B_cell,Endothelial_Cell,Fibroblast,Macrophage,Monocyte_Cell,MP,Muscle_Cell,myoFB,NK_Cell,Smooth_Muscle_Cell,T_cell,unknown)
new.cluster.ids <- c(rep("AP",length(Adipocyte_Adipose_Tissue)),
                     rep("B",length(B_cell)),
                     rep("EC",length(Endothelial_Cell)),
                     rep("FB",length(Fibroblast)),
                     rep("mDC",length(Macrophage)),
                     rep("Mono",length(Monocyte_Cell)),
                     rep("MP",length(MP)),
                     rep("muscle",length(Muscle_Cell)),
                     rep("myoFB",length(myoFB)),
                     rep("NK",length(NK_Cell)),
                     rep("SMC",length(Smooth_Muscle_Cell)), 
                     rep("T",length(T_cell)),
                     rep("Unknown",length(unknown))
                     )
scRNA@meta.data$cellType <- plyr::mapvalues(x = scRNA$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
# scRNA <- subset(scRNA, pre_cellType %in% c("T cell","B cell","Myeloid cell","Epithelial cell","Fibrocyte","Endothelial cell","Mast cell",'MIX'))

pdf(file="p_cellType_anno2.pdf",height = 10, width = 7)
DimPlot(scRNA,group.by = "cellType",label=T)
dev.off()

scRNA@meta.data$cellType=factor(scRNA@meta.data$cellType,levels =c("AP","B","mDC","EC","FB","MP","Mono","muscle","NK","myoFB","SMC","T","Unknown"))
pdf(file="p_genedotplot-celltype2.pdf",height = 15, width = 15)
DotPlot(scRNA,group.by = 'cellType', features = unique(genes_to_check)) + coord_flip()
dev.off()
# 提取个样本数据画dimplot
scRNASe=subset(x = scRNA, subset = (orig.ident == "Se"))
scRNAJe=subset(x = scRNA, subset = (orig.ident == "Je"))
scRNASi=subset(x = scRNA, subset = (orig.ident == "Si"))
scRNAJi=subset(x = scRNA, subset = (orig.ident == "Ji"))
scRNAST=subset(x = scRNA, subset = (orig.ident == "ST"))
scRNAJT=subset(x = scRNA, subset = (orig.ident == "JT"))
p2=DimPlot(scRNASe,group.by = "cellType",label=T)
ggsave(filename = "p_DimPlot_Se.pdf", plot = p2, height = 10, width = 7)
p2=DimPlot(scRNAJe,group.by = "cellType",label=T)
ggsave(filename = "p_DimPlot_Je.pdf", plot = p2, height = 10, width = 7)
p2=DimPlot(scRNASi,group.by = "cellType",label=T)
ggsave(filename = "p_DimPlot_Si.pdf", plot = p2, height = 10, width = 7)
p2=DimPlot(scRNAJi,group.by = "cellType",label=T)
ggsave(filename = "p_DimPlot_Ji.pdf", plot = p2, height = 10, width = 7)
p2=DimPlot(scRNAST,group.by = "cellType",label=T)
ggsave(filename = "p_DimPlot_ST.pdf", plot = p2, height = 10, width = 7)
p2=DimPlot(scRNAJT,group.by = "cellType",label=T)
ggsave(filename = "p_DimPlot_JT.pdf", plot = p2, height = 10, width = 7)

####各细胞类型占比####
library(Seurat)
library(ggplot2)
library(dplyr)
# scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c("Se","Je","Si","Ji","ST","JT"))
# scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c("Se","Si","ST","Je","Ji","JT"))
scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c("Je","Se","Ji","Si","JT","ST"))
Ratio <- scRNA@meta.data %>%group_by(orig.ident,cellType) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(Freq = n/sum(n)*100)

p1=ggplot(Ratio, aes(x = orig.ident, y = Freq, fill = cellType))+
  geom_col()+
  geom_text(aes(label = paste(round(Freq, 1),"%")), 
            position = position_stack(vjust = 0.5))+
  theme_classic()+
  scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                              "#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                              "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02"))+ theme(text = element_text(size = 18))+theme(legend.text = element_text(size = 18),axis.text.x= element_text(size = 18),axis.text.y= element_text(size = 18))
ggsave(filename = "p_Ratio2.pdf", plot = p1, height = 10, width = 7)
Ratio <- scRNA@meta.data %>%group_by(orig.ident,cellType) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(Freq = n/sum(n)*100)
write.csv(file='ratio.csv',Ratio)
p1=ggplot(Ratio, aes(x = orig.ident, y = Freq, fill = cellType))+
  geom_col()+
  geom_text(aes(label = paste(round(Freq, 1),"%")), 
            position = position_stack(vjust = 0.5))+
  theme_classic()+
  scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                              "#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                              "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02"))
ggsave(filename = "p_Ratio.pdf", plot = p1, height = 10, width = 7)
# scRNA<-readRDS("JS_seurat_cluster_filt_Anchors.rds")
features<-c('Adipoq','Pdgfra','Msln','Jam2','Prox1','Steap4','Myocd','Mafb','Cybb','Flt3','Cpa3','Csf3r','Ms4a1','Klrd1','Il7r','Dcdc2a','Erbb4')
pdf(file="p_FeaturePlot_Anchors.pdf",height = 20, width = 30)
FeaturePlot(object = scRNA, features = features)
dev.off()

# 不加label
scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c("Se","Je","Si","Ji","ST","JT"))
# scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c("Se","Si","ST","Je","Ji","JT"))
Ratio <- scRNA@meta.data %>%group_by(orig.ident,cellType) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(Freq = n/sum(n)*100)
write.csv(file='ratio1.csv',Ratio)
p1=ggplot(Ratio, aes(x = orig.ident, y = Freq, fill = cellType))+
  geom_col()+
  theme_classic()+
  scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                              "#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                              "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02"))
ggsave(filename = "p_Ratio2.pdf", plot = p1, height = 5, width = 3.5)


# find all markers of cluster 2
cluster9.markers <- FindMarkers(scRNA, ident.1 = 9,only.pos = TRUE)
cluster15.markers <- FindMarkers(scRNA, ident.1 = 15,only.pos = TRUE)
cluster16.markers <- FindMarkers(scRNA, ident.1 = 16,only.pos = TRUE)
cluster20.markers <- FindMarkers(scRNA, ident.1 = 20,only.pos = TRUE)
cluster30.markers <- FindMarkers(scRNA, ident.1 = 30,only.pos = TRUE)
cluster33.markers <- FindMarkers(scRNA, ident.1 = 33,only.pos = TRUE)
cluster36.markers <- FindMarkers(scRNA, ident.1 = 36,only.pos = TRUE)

cluster12.markers <- FindMarkers(scRNA, ident.1 = 12,only.pos = TRUE)

cluster10.markers <- FindMarkers(scRNA, ident.1 = 10,only.pos = TRUE)
cluster11.markers <- FindMarkers(scRNA, ident.1 = 11,only.pos = TRUE)
cluster22.markers <- FindMarkers(scRNA, ident.1 = 22,only.pos = TRUE)
cluster28.markers <- FindMarkers(scRNA, ident.1 = 28,only.pos = TRUE)


head(cluster9.markers, n = 10)
head(cluster15.markers, n = 10)
head(cluster16.markers, n = 10)
head(cluster20.markers, n = 10)
head(cluster30.markers, n = 10)
head(cluster33.markers, n = 10)
head(cluster36.markers, n = 10)

head(cluster12.markers, n = 10)

head(cluster10.markers, n = 10)
head(cluster11.markers, n = 10)
head(cluster22.markers, n = 10)
head(cluster28.markers, n = 10)




# 展示细胞簇的基因表达情况

genes_to_check = c('Adipoq', 'Cfd','Gsn','Camk1d','Hp',#adpocyte markers
                   'Pdgfra', 'Hmcn1','Col8a1','Col6a5','Slit2',#ASPC markers
                   'Msln', 'Nkain2','Nav2','Trdn','Ppfibp2',#mesothelial markers
                   'Mgp','Ptprj','Cytl1','Eln','Pcdh9',#vascular markers
                   'Mreg','Il12b','Strip2','Cacnb3','Traf1',#immune markers
                   'Prlr','Nrxn3','Agps','Tmem56','Ptn',#epithelial markers
                   'Jam2',#Endothelial
                   'Prox1',#Lymphatic Endo
                   'Steap4',#Pericyte
                   'Myocd',#Smooth Muscle
                   'Mafb',#Macrophage
                   'Cybb',#Monocyte
                   'Flt3',#Dendritic Cell
                   'Cpa3',#Mast Cell
                   'Csf3r',#Neutrophil
                   'Ms4a1',#B Cell
                   'Klrd1',#NK Cell
                   'Il7r',#T Cell
                   'Dcdc2a',#Male Epithelial
                   'Erbb4'#Female Epithelial
                )

genes_to_check = c('Adipoq','Cfd','Gsn','Camk1d','Hp',#adpocyte 5,29
                   'Pdgfra',#ASPC   6,14
                   'Msln',#mesothelium
                   'Jam2',#Endothelial  4,24,31
                   'Prox1',#Lymphatic Endo
                   'Steap4',#Pericyte   17,27
                   'Myocd',#Smooth Muscle,21
                   'Mafb',#Macrophage   19
                   'Cybb',#Monocyte 3
                   'Flt3',#Dendritic Cell   34
                   'Cpa3',#Mast Cell    
                   'Csf3r',#Neutrophil
                   'Ms4a1',#B Cell  10
                   'Klrd1',#NK Cell 22
                   'Il7r',#T Cell   33,11,28
                   'Dcdc2a',#Male Epithelial    
                   'Erbb4'#Female Epithelial    32
                )

genes_to_check = c('Adipoq','Cfd','Gsn','Camk1d','Hp'#adpocyte
                )

genes_to_check = c('Adipoq', 'Cfd','Gsn','Camk1d','Hp',#adpocyte markers
                   'Pdgfra', 'Hmcn1','Col8a1','Col6a5','Slit2',#ASPC markers
                   'Msln', 'Nkain2','Nav2','Trdn','Ppfibp2',#mesothelial markers
                   'Mgp','Ptprj','Cytl1','Eln','Pcdh9',#vascular markers
                   'Mreg','Il12b','Strip2','Cacnb3','Traf1',#immune markers
                   'Prlr','Nrxn3','Agps','Tmem56','Ptn'#epithelial markers
                )

pdf(file="p_clustergenedotplot-adpocyte.pdf",height = 15, width = 15)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(genes_to_check),cluster.idents = T) + coord_flip()
dev.off()


####method-Je/Se####
library(Seurat)
library(tidyverse)
library(patchwork)
# 导入单样本文件
dir = c(
  "./Je/filtered_feature_bc_matrix",
  "./Se/filtered_feature_bc_matrix"
)

# 读取并创建Seurat Object
scRNAlist <- list()  # 创建一个空的列表,其中包含三个seurat对象
names(dir) = c("Je", 'Se')
sample_name = c("Je", 'Se')
for (i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = sample_name[i])
}

# 合并Seurat object
# scRNA <- merge(scRNAlist[[1]],y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]]))
# dim(scRNA)

#  数据标准化和选择高变基因
# 每一个样本分别进行数据标准化和提取高变基因
features <- SelectIntegrationFeatures(object.list = scRNAlist)
for (i in 1:length(scRNAlist)){
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-")
    scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
    scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
    scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]],selection.method = 'vst',nfeatures = 3000) # 这里我们只需要2000的基因
    scRNAlist[[i]] <- ScaleData(scRNAlist[[i]], features = features, verbose = FALSE)
    scRNAlist[[i]] <- RunPCA(scRNAlist[[i]], features = features, verbose = FALSE)
}

scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, dims = 1:30)
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)

DefaultAssay(scRNA) <- "integrated"

# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, verbose = FALSE)
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
# pdf(file="p_ElbowPlot_Anchors.pdf",,height = 7, width = 12)
# ElbowPlot(scRNA,ndims = 30)
# dev.off()
# t-SNE and Clustering
# scRNA <- RunTSNE(scRNA, reduction = "pca", dims = 1:30)
scRNA  <- RunUMAP(scRNA , dims=1:30,reduction='pca')
scRNA <- FindNeighbors(scRNA, reduction = "pca", dims = 1:30)
scRNA <- FindClusters(scRNA,resolution = 0.8)
# saveRDS(scRNA,"scRNA_rmMT_update.seurat.RDS")
# 查看各样本信息
table(scRNA@meta.data$orig.ident)
p1=DimPlot(scRNA,group.by = "orig.ident")
# 查看聚类信息
# p2=DimPlot(scRNA,group.by = "seurat_clusters",label=T)
# combine<-CombinePlots(list(p1,p2))
ggsave(filename = "p_DimPlot_sample-JESE.pdf", plot = p1, height = 10, width = 14)

####method-N7####
library(Seurat)
library(tidyverse)
library(patchwork)
# 导入单样本文件
dir = c(
  "./Ji/filtered_feature_bc_matrix",
  "./Si/filtered_feature_bc_matrix"
)

# 读取并创建Seurat Object
scRNAlist <- list()  # 创建一个空的列表,其中包含三个seurat对象
names(dir) = c( "Ji", 'Si')
sample_name = c( "Ji", 'Si')
for (i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = sample_name[i])
}

# 合并Seurat object
# scRNA <- merge(scRNAlist[[1]],y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]]))
# dim(scRNA)

#  数据标准化和选择高变基因
# 每一个样本分别进行数据标准化和提取高变基因
features <- SelectIntegrationFeatures(object.list = scRNAlist)
for (i in 1:length(scRNAlist)){
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-")
    scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
    scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
    scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]],selection.method = 'vst',nfeatures = 3000) # 这里我们只需要2000的基因
    scRNAlist[[i]] <- ScaleData(scRNAlist[[i]], features = features, verbose = FALSE)
    scRNAlist[[i]] <- RunPCA(scRNAlist[[i]], features = features, verbose = FALSE)
}

scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, dims = 1:30)
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)

DefaultAssay(scRNA) <- "integrated"

# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, verbose = FALSE)
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
# pdf(file="p_ElbowPlot_Anchors.pdf",,height = 7, width = 12)
# ElbowPlot(scRNA,ndims = 30)
# dev.off()
# t-SNE and Clustering
# scRNA <- RunTSNE(scRNA, reduction = "pca", dims = 1:30)
scRNA  <- RunUMAP(scRNA , dims=1:30,reduction='pca')
scRNA <- FindNeighbors(scRNA, reduction = "pca", dims = 1:30)
scRNA <- FindClusters(scRNA,resolution = 0.8)
# saveRDS(scRNA,"scRNA_rmMT_update.seurat.RDS")
# 查看各样本信息
table(scRNA@meta.data$orig.ident)
p1=DimPlot(scRNA,group.by = "orig.ident")
# 查看聚类信息

# p2=DimPlot(scRNA,group.by = "seurat_clusters",label=T)
# combine<-CombinePlots(list(p1,p2))
ggsave(filename = "p_DimPlot_sample_JISI.pdf", plot = p1, height = 10, width = 14)

######JTST#####
library(Seurat)
library(tidyverse)
library(patchwork)
# 导入单样本文件
dir = c(
  "./JT/filtered_feature_bc_matrix",
  "./ST/filtered_feature_bc_matrix"
)

# 读取并创建Seurat Object
scRNAlist <- list()  # 创建一个空的列表,其中包含三个seurat对象
names(dir) = c("JT", 'ST')
sample_name = c("JT", 'ST')
for (i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = sample_name[i])
}

# 合并Seurat object
# scRNA <- merge(scRNAlist[[1]],y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]]))
# dim(scRNA)

#  数据标准化和选择高变基因
# 每一个样本分别进行数据标准化和提取高变基因
features <- SelectIntegrationFeatures(object.list = scRNAlist)
for (i in 1:length(scRNAlist)){
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-")
    scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
    scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
    scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]],selection.method = 'vst',nfeatures = 3000) # 这里我们只需要2000的基因
    scRNAlist[[i]] <- ScaleData(scRNAlist[[i]], features = features, verbose = FALSE)
    scRNAlist[[i]] <- RunPCA(scRNAlist[[i]], features = features, verbose = FALSE)
}

scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, dims = 1:30)
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)

DefaultAssay(scRNA) <- "integrated"

# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, verbose = FALSE)
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
# pdf(file="p_ElbowPlot_Anchors.pdf",,height = 7, width = 12)
# ElbowPlot(scRNA,ndims = 30)
# dev.off()
# t-SNE and Clustering
# scRNA <- RunTSNE(scRNA, reduction = "pca", dims = 1:30)
scRNA  <- RunUMAP(scRNA , dims=1:30,reduction='pca')
scRNA <- FindNeighbors(scRNA, reduction = "pca", dims = 1:30)
scRNA <- FindClusters(scRNA,resolution = 0.8)
# saveRDS(scRNA,"scRNA_rmMT_update.seurat.RDS")
# 查看各样本信息
table(scRNA@meta.data$orig.ident)
p1=DimPlot(scRNA,group.by = "orig.ident")
# 查看聚类信息
# p2=DimPlot(scRNA,group.by = "seurat_clusters",label=T)
# combine<-CombinePlots(list(p1,p2))
ggsave(filename = "p_DimPlot_sample-JTST.pdf", plot = p1, height = 10, width = 14)


features<-c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_FeaturePlot_Adipocyte.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_Adipocyte.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()

features<-c('Cd68','Fcgr1','Naaa','Lyz2','Ccl12','Sepp1')
pdf(file="p_FeaturePlot_Macrophage.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_Macrophage.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()

features<-c('Cd93','Vwf','Emcn','Egfl7','Flt1','Id3','Vwf','Ace','Cd62p','Adamts13','Pecam1','Vcam1','Icam1','Icam2','Cd47','Sele','Cdh5','Nectin2','Esam','Itgb3','Cd151','Cd248','Mcam','Itga4')
pdf(file="p_FeaturePlot_Endothelial.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_Endothelial.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()

features<-c('Gfap','Slc1a2','Acsl6','Agt','Aqp4','Apoe','S100b','Sox9')
pdf(file="p_FeaturePlot_Astrocyte.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_Astrocyte.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()

features<-c('Pxk','Ms4a1','Cd19','Cd74','Cd79a','Ighd')
pdf(file="p_FeaturePlot_B_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_B_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()

features<-c('Vim','Pdgfrb','Lum','Col6a2','Vtn','Mfap5','Cd13','Cd90','Dlk1','Cd26','Fap','Te-7','S100a4','1b10','Cd29','Mas516','Fdgfr','Hsp47','P4hb')
pdf(file="p_FeaturePlot_Fibroblast.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_Fibroblast.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()

# features<-c('Itgam','Cx3cr1','Tmem119','P2ry12','Aif1','Csf1r')
features<-c('lba1','Cx3cr1','Tmem119','P2ry12','Cd11b','Csf1r')
pdf(file="p_FeaturePlot_Microglial.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_Microglial.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()

features<-c('Itgax','Zbtb46','Cd86','Lamp3','Cd83','Cd1a')
pdf(file="p_FeaturePlot_DendriticCell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_DendriticCell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()


features<-c('Trdc','Nkg7','Klrf1','Klrd1','Gnly','Ncr1')
pdf(file="p_FeaturePlot_NK_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_NK_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()


features<-c('Nes','Neurod1','Eomes')
pdf(file="p_FeaturePlot_Neural_stem_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_Neural_stem_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()


features<-c('Acta2','Myl9','Rgs5','Mylk','Hhip','Myh11')
pdf(file="p_FeaturePlot_SmoothMuscle_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_SmoothMuscle_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()


features<-c('Trbc2','Cd3d','Cd3g','Cd3e','Il7r','Ltb')
pdf(file="p_FeaturePlot_T_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_T_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()



features<-c('Tnnt3','Ttn','Myod1','Actc1','Mylpf','Tnnt1','Actal')
pdf(file="p_FeaturePlot_Muscle_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_Muscle_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()


features<-c('A2b5','Aspa','Bmp4','Cd82','Cd9','Cnp','Cspg4','Enpp4','Mag','Mal','Mog','Nkx2-2','Nkx6-2','Olig1','Olig2','Pdgfra','Plp1','Smarca4','Sox10','Sox17','Tfr','Tmem10','Zfp191','Zfp488','Zfp536')
pdf(file="p_FeaturePlot_Oligodendrocyte_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_Oligodendrocyte_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()

features<-c('Abhd5','Cd49f','Cnn1','Ctp1a','Ctp2','Connexins','Eaat1','Eaat2','Fbp','Fmn2','G6pc','Gat-1','Gat-3','Gfap','Gk','Gp','Gpat3','Gs','Klr4.1','Ldh5','Mct4','Nebl','Nf1a','Nf1b','Nf1x','Nkx2-1','Nkx3-1','Nkx6-1','Pdk4','Pdlim7','Pdz','Pfk','Pfkfb3','Pkm2','S100b','Sox9','Synpo2')
pdf(file="p_FeaturePlot_Astrocyte_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_Astrocyte_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()

features<-c('B220','Btla','Cadm1','Cd14','Cd163','Cd1a','Cd1c','Cd2','Cd207','Cd209','Cd45ra','Cd83','Cd8a','Clec10a','Clec4c','Clec9a','Cx3cr1','Fcer1','Id2','Ifn-α','Itgae','Itgam','Itgax','Lilrb4','Ly75','Mrc1','Notch2','Nrp1','Sirpa','Siglech','Thbd','Xcr1')
pdf(file="p_FeaturePlot_Dendritic_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_Dendritic_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()


features<-c('E2f7','Top2a','Cdca2','Ckap2l','Kif23','Kif11','Hmmr','Mki67','Smc4','Cenpa','Ndc80','Kif15','Cit','Aspm','Anln','Cenpe','Cenpf','Diaph3','Prc1','Nusap1','Bub1','Cep55')
pdf(file="p_FeaturePlot_VascularStem_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_VascularStem_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()



features<-c('')
pdf(file="p_FeaturePlot_B_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_B_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()


features<-c('')
pdf(file="p_FeaturePlot_B_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_B_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()



features<-c('')
pdf(file="p_FeaturePlot_B_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_B_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()

features<-c('')
pdf(file="p_FeaturePlot_B_cell.pdf",height = 50, width = 25)
FeaturePlot(object = scRNA, features = features)
dev.off()
# genes_to_check = c('Adipoq','Pnpla2','Plin1','Cidec','Apoc1','Fabp4')
pdf(file="p_dotplot_B_cell.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(features),cluster.idents = T) + coord_flip()
dev.off()

genes_to_check = c("Cidec","Ghr","Mmd","Nnat","Sgcd","A530053G22Rik","Acsl1","Car3","Ctcflos","Faml3a","Mpp7","Npr3","Nrg4","Prkd1","Slc1a3","Slc7a10","Sntb1","Sorbs1","Tenm4","Aff3","Arhgap15","Bach2","Bcl1la","Blnk","Dock2","Elmo1","Gm26740","Ikzf1","Ikzf3","Inpp4b","Lyn","Mef2c","Mndal","Prkcb","Ptprc","Ralgps2","Ripor2","Slc9a7","St8sia4","Syk","Trim30a","Parp8","Plxnc1","Samsnl","Tmem150c","Tmtc2","Ankrd33b","Cacnb3","Cd83","Ly75","Mreg","Traf1","Ptprb","Adgrf5","Aqp1","Cyyr1","Flt1","Sema6a","Cd200","Etl4","Gpihbp1","Ldb2","Mecom","Pdgfd","Pecam1","Pkp4","Adamtsl1","Dcn","Ddr2","Ebf2","Egfr","Fstl1","Naaladl2","Tmeff2","Tnxb","Adgrd1","Antxr1","Bicc1","Clec3b","Col1a1","Col1a2","Col3a1","Col5a2","Flrt2","Cd74","H2-Aa","H2-Ab1","H2-Eb1","Mctp1","Abca9","Adgrel","Ctsc","Ehd4","F13a1","Ly86","Mrc1","Nrros","Ptprj","Tpm1","Acta1","Actn3","Ank3","Ano5","Atp2a1","Ckm","Cmya5","Fgf13","Kcnma1","Ldb3","Mylk4","Myom1","Neb","0bscn","Pde4dip","Pvalb","Pygm","Ryr1","Tnni2","Tnnt3","Trdn","Ttn","Adgrl3","Ank3","Ano4","Cadm2","Lrp1b","Map2","Ncaml","Nrxnl","Specc1","Trdc","Cacna1c","Cald1","Gucy1a1","Gucy1a2","Myh11","Notch3","Pde3a","Prkg1","Sgip1","Stac","Timp3","Trbc2","Cd34","Cd31","Cd133","Vegfr2","Vwf","Uea-1","Sm-mhc","Sca-1","Vegfr1","Sma","Pdgfr-b","desmin","endoglin","Ng2","Notch3","Rgs5","Cd146","Cd4","Cd73","Cd90","Cd105","Cd166","Stro-1","Ssea1")
pdf(file="p_dotplot_allImarkers.pdf",,height = 30, width = 30)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(genes_to_check),cluster.idents = T) + coord_flip()
dev.off()


scRNAAP=subset(x = scRNA, subset = (cellType == "AP"))
scRNASeJeAP=subset(x = scRNAAP, subset = (orig.ident == "Se"|orig.ident == "Je"))
table(Idents(scRNASeJeAP))

pdf(file="p_DimPlot_AP_SeJe.pdf",height = 20, width = 10)
DimPlot(scRNASeJeAP, reduction = "umap",group.by = "seurat_clusters",label=T,label.size = 16)+ xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = 40), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))+theme(axis.text.x=element_text(vjust=1,size=40))+theme(axis.text.y=element_text(vjust=1,size=40))+theme(plot.title=element_text(size=40))
dev.off()

Idents(scRNASeJeAP)='seurat_clusters'
DefaultAssay(scRNASeJeAP) <- "RNA"

scRNASeJeAP.markers <- FindAllMarkers(scRNASeJeAP, only.pos = TRUE)
write.csv(file='scRNASeJeAP.markers.csv',scRNASeJeAP.markers)

pdf(file="p_VlnPlot_AP_brownall_SeJe.pdf",height = 10, width = 15)
VlnPlot(scRNASeJeAP, features = unique(c("Col5a3","Cxcl14","Bmper","DPP4","Pi16","Clec11a","Gdf10","Trpv1","Sca1","Cd81","Vipr2","Gli1","Pparg","Trpv1","Myh11","Adipoq","Pparg","Vegfa","Ly6a","Cd34","Thy1","Mfap5","Clec3b","Ecm1","MMP12","MMP19","Fn1","Mrc1","Clec10a","C1qa","C1qb","CD36","Lpl","Lipa","Fabp4","Trem2","Plin2","Treml4","Plac8","Ly6c2")),ncol=6)
dev.off()