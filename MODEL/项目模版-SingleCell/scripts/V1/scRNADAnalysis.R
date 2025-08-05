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
p_QC <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "p_QC.pdf", height = 6, width = 6, plot = p_QC)
# 过滤
scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
# 不去批次效应
# Normalization
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA,nfeatures = 3000)
scRNA <- ScaleData(scRNA,vars.to.regress = 'percent.mt')
scRNA <- RunPCA(scRNA)
scRNA <- FindNeighbors(scRNA, reduction = "pca", dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:30) 
# 查看各样本信息
table(scRNA@meta.data$orig.ident)
p1=DimPlot(scRNA,group.by = "orig.ident")
# 查看聚类信息
p2=DimPlot(scRNA,group.by = "seurat_clusters",label=T)
combine<-CombinePlots(list(p1,p2))
ggsave(filename = "p_DimPlotNoharmony.pdf", plot = combine, height = 10, width = 14)



# harmony去除批次效应
scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident")
scRNA  <- FindNeighbors(scRNA , dims = 1:30 , reduction = "harmony")
scRNA  <- FindClusters(scRNA, save.snn=T , resolution = 0.8)
scRNA  <- RunUMAP(scRNA , dims=1:30,reduction='harmony')
table(scRNA@meta.data$orig.ident)
p1=DimPlot(scRNA,group.by = "orig.ident")
# 查看聚类信息
p2=DimPlot(scRNA,group.by = "seurat_clusters",label=T)
combine<-CombinePlots(list(p1,p2))
ggsave(filename = "p_DimPlotharmony.pdf", plot = combine, height = 10, width = 14)





# Anchors去除批次效应
#  数据标准化和选择高变基因
# 每一个样本分别进行数据标准化和提取高变基因
features <- SelectIntegrationFeatures(object.list = scRNAlist)
for (i in 1:length(scRNAlist)){
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-")
    scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
    scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
    scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]],selection.method = 'vst',nfeatures = 3000) # 这里我们只需要3000的基因
    scRNAlist[[i]] <- ScaleData(scRNAlist[[i]], features = features, verbose = FALSE)
    scRNAlist[[i]] <- RunPCA(scRNAlist[[i]], features = features, verbose = FALSE)
}

# 寻找Integration anchors
# 我们使用PCA来寻找anchor
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, dims = 1:30)
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)

DefaultAssay(scRNA) <- "integrated"

# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, verbose = FALSE)
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
pdf(file="p_ElbowPlotAnchors.pdf",height = 6, width = 6)
ElbowPlot(scRNA,ndims = 30)
dev.off()
# t-SNE and Clustering
# scRNA <- RunTSNE(scRNA, reduction = "pca", dims = 1:30)
scRNA  <- RunUMAP(scRNA , dims=1:30,reduction='pca')
scRNA <- FindNeighbors(scRNA, reduction = "pca", dims = 1:30)
scRNA <- FindClusters(scRNA,resolution = 0.8)
# saveRDS(scRNA,"scRNA.seurat.RDS")

# 查看各样本信息
table(scRNA@meta.data$orig.ident)
p1=DimPlot(scRNA,group.by = "orig.ident")
# 查看聚类信息
p2=DimPlot(scRNA,group.by = "seurat_clusters",label=T)
combine<-CombinePlots(list(p1,p2))
ggsave(filename = "p_DimPlotAbchors.pdf", plot = combine, height = 10, width = 14)

DefaultAssay(scRNA) <- "RNA"
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
pdf(file="p_genedotplot-celltype2.pdf",height = 10, width = 7)
DotPlot(scRNA,group.by = 'cellType', features = unique(genes_to_check)) + coord_flip() +RotatedAxis()+ theme(text = element_text(size = 16))
dev.off()





