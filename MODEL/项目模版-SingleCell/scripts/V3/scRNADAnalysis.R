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

saveRDS(file='scRNA_rawdata.rds',scRNA)
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
ggsave(filename = "p_DimPlotNoharmony1.pdf", plot = p1, height = 6, width = 6)
# 查看聚类信息
p2=DimPlot(scRNA,group.by = "seurat_clusters",label=T)
ggsave(filename = "p_DimPlotNoharmony2.pdf", plot = p2, height = 6, width = 6)
combine<-CombinePlots(list(p1,p2))
ggsave(filename = "p_DimPlotNoharmony.pdf", plot = combine, height = 6, width = 12)
saveRDS(file='scRNA_filt.unharmony.rds',scRNA)


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
ggsave(filename = "p_DimPlotharmony.pdf", plot = combine, height = 6, width = 6)
saveRDS(file='scRNA_filt.harmony.rds',scRNA)





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
saveRDS(scRNA,"scRNA.filt.anchors.rds")

# 查看各样本信息
table(scRNA@meta.data$orig.ident)
scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c("Je","Ji","JT","Se","Si","ST"))
p1=DimPlot(scRNA,group.by = "orig.ident")
ggsave(filename = "p_DimPlotAnchors1.pdf", plot = p1, height = 6, width = 6)
# 查看聚类信息

p2=DimPlot(scRNA,group.by = "seurat_clusters",label=T)
ggsave(filename = "p_DimPlotAnchors2.pdf", plot = p2, height = 6, width = 6)
combine<-CombinePlots(list(p1,p2))
ggsave(filename = "p_DimPlotAnchors.pdf", plot = combine, height = 6, width = 12)

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
saveRDS(scRNA,"scRNA.celltype.rds")
scRNA=readRDS('scRNA_seurat_cellType.RDS')

pdf(file="p_cellType_anno.pdf",height = 6, width = 6)
DimPlot(scRNA,group.by = "cellType",label=T)
dev.off()

pdf(file="p_cellType_anno.split.pdf",height = 6, width = 36)
DimPlot(scRNA,group.by = "cellType",label=T,split.by = "orig.ident")
dev.off()

# 提取个样本数据画dimplot
scRNASe=subset(x = scRNA, subset = (orig.ident == "Se"))
scRNAJe=subset(x = scRNA, subset = (orig.ident == "Je"))
scRNASi=subset(x = scRNA, subset = (orig.ident == "Si"))
scRNAJi=subset(x = scRNA, subset = (orig.ident == "Ji"))
scRNAST=subset(x = scRNA, subset = (orig.ident == "ST"))
scRNAJT=subset(x = scRNA, subset = (orig.ident == "JT"))
p2=DimPlot(scRNASe,group.by = "cellType",label=T)
ggsave(filename = "p_DimPlot_Se.pdf", plot = p2, height = 6, width = 6)
p2=DimPlot(scRNAJe,group.by = "cellType",label=T)
ggsave(filename = "p_DimPlot_Je.pdf", plot = p2, height = 6, width = 6)
p2=DimPlot(scRNASi,group.by = "cellType",label=T)
ggsave(filename = "p_DimPlot_Si.pdf", plot = p2, height = 6, width = 6)
p2=DimPlot(scRNAJi,group.by = "cellType",label=T)
ggsave(filename = "p_DimPlot_Ji.pdf", plot = p2, height = 6, width = 6)
p2=DimPlot(scRNAST,group.by = "cellType",label=T)
ggsave(filename = "p_DimPlot_ST.pdf", plot = p2, height = 6, width = 6)
p2=DimPlot(scRNAJT,group.by = "cellType",label=T)
ggsave(filename = "p_DimPlot_JT.pdf", plot = p2, height = 6, width = 6)


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

scRNA@meta.data$cellType=factor(scRNA@meta.data$cellType,levels =c("AP","B","mDC","EC","FB","MP","Mono","muscle","NK","myoFB","SMC","T","Unknown"))
pdf(file="anno_p_genedotplot-cellType.pdf",height = 12, width = 6)
DotPlot(scRNA,group.by = 'cellType', features = unique(genes_to_check)) + coord_flip() +RotatedAxis()+ theme(text = element_text(size = 16))+theme(axis.text.x = element_text(size=20))+theme(axis.text.y = element_text(size=20))
dev.off()

scRNA@meta.data$seurat_clusters=factor(scRNA@meta.data$seurat_clusters,levels =c("0","1","2","7","8","27","10","28","4","5","24","26","31","6","13","14","25","3","18","19","33","34","12","29","22","17","21","9","15","16","20","11","30","32","35","36","23"))
pdf(file="p_clustersgenedotplot2.pdf",height = 6, width = 6)
DotPlot(scRNA,group.by = 'seurat_clusters', features = unique(genes_to_check)) + coord_flip()
dev.off()

####各细胞类型占比####
library(ggplot2)
library(dplyr)
# scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c("Se","Je","Si","Ji","ST","JT"))
# scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c("Se","Si","ST","Je","Ji","JT"))
scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c("Je","Se","Ji","Si","JT","ST"))
Ratio <- scRNA@meta.data %>%group_by(orig.ident,cellType) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(Freq = n/sum(n)*100)

# 添加label
p1=ggplot(Ratio, aes(x = orig.ident, y = Freq, fill = cellType))+
  geom_col()+
  geom_text(aes(label = paste(round(Freq, 1),"%")), 
            position = position_stack(vjust = 0.5))+
  theme_classic()+
  scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                              "#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                              "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02"))
ggsave(filename = "p_Ratio.pdf", plot = p1, height = 6, width = 6)

# 不添加label
p1=ggplot(Ratio, aes(x = orig.ident, y = Freq, fill = cellType))+
  geom_col()+
  theme_classic()+
  scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                              "#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                              "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02"))
ggsave(filename = "p_Ratio.pdf", plot = p1, height = 6, width = 6)


#####按样本拆分scRNA#####
scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))
scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
# scRNASTJT=subset(x = scRNA, subset = (orig.ident == "ST"|orig.ident == "JT"))

pdf(file="Nucb2_scRNASeJe.pdf",height = 6, width = 12)
FeaturePlot(object = scRNASeJe, features = "Nucb2",split.by = 'orig.ident')
dev.off()

pdf(file="Nucb2_scRNASiJi.pdf",height = 6, width = 12)
FeaturePlot(object = scRNASiJi, features = "Nucb2",split.by = 'orig.ident')
dev.off()

pdf(file="Nucb2_scRNASTJT.pdf",height = 6, width =12)
FeaturePlot(object = scRNASTJT, features = "Nucb2",split.by = 'orig.ident')
dev.off()





#合并两样本dimplot,展示gene
p1=DimPlot(scRNASeJe,group.by = 'orig.ident',
           cols = c("grey","grey"),
           pt.size = 1)+NoLegend()+ labs(title="Nucb2")

gene1 <- "Nucb2"

pos_gene=scRNASeJe@reductions$umap@cell.embeddings

pos_gene1=pos_gene[scRNASeJe@assays$RNA@data[gene1,]>0,]

pos_gene1Se=pos_gene1[grep('Se',rownames(pos_gene1)),]
pos_gene1Je=pos_gene1[grep('Je',rownames(pos_gene1)),]

p=p1+geom_point(data = as.data.frame(pos_gene1Se),
              aes(x=UMAP_1,y=UMAP_2), 
              shape = 16,size = 1.4, color="#6a5acd")+
  geom_point(data = as.data.frame(pos_gene1Je),
             aes(x=UMAP_1,y=UMAP_2),
             shape = 16,size = 1.4, color="#ff69b4")+theme(text = element_text(size = 30))+theme(axis.text.x = element_text(size=20),
           axis.text.y = element_text(size=20))
ggsave(filename = "SeJe-Dimplotgene.pdf", height = 6, width = 6, plot =p)

#合并两样本dimplot,展示gene
p1=DimPlot(scRNASiJi,group.by = 'orig.ident',
           cols = c("grey","grey"),
           pt.size = 1)+NoLegend()+ labs(title="Nucb2")

gene1 <- "Nucb2"

pos_gene=scRNASiJi@reductions$umap@cell.embeddings

pos_gene1=pos_gene[scRNASiJi@assays$RNA@data[gene1,]>0,]

pos_gene1Se=pos_gene1[grep('Si',rownames(pos_gene1)),]
pos_gene1Je=pos_gene1[grep('Ji',rownames(pos_gene1)),]

p=p1+geom_point(data = as.data.frame(pos_gene1Se),
              aes(x=UMAP_1,y=UMAP_2), 
              shape = 16,size = 1.4, color="#6a5acd")+
  geom_point(data = as.data.frame(pos_gene1Je),
             aes(x=UMAP_1,y=UMAP_2),
             shape = 16,size = 1.4, color="#ff69b4")+theme(text = element_text(size = 30))+theme(axis.text.x = element_text(size=20),
           axis.text.y = element_text(size=20))
ggsave(filename = "SiJi-Dimplotgene.pdf", height = 6, width = 6, plot =p)



# 通路活性分析
library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(SeuratData)

fatty <- read.gmt("fat.gmt")
head(fatty$term)
head(fatty$gene)
geneSets <- lapply(unique(fatty$term), function(x){print(x);fatty$gene[fatty$term == x]})
names(geneSets) <- unique(fatty$term)
features=rownames(unique(fatty$term))
features=gsub( "_", "-",features)#后续要加入到seurat，行名不能含有下划线

obj=scRNASeJe # 根据不同的对象进行调整
sample1=''
sample2=''

cells_rankings <- AUCell_buildRankings(obj@assays$RNA@data, plotStats=TRUE) 
# cells_rankings
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
aucs <- getAUC(cells_AUC)
# 将aucs加入seurat
# aucs=t(aucs)
obj_AUC <- CreateAssayObject(counts = aucs)
obj@assays$AUC <- obj_AUC
obj@assays$AUC@key <- "rna_"
DefaultAssay(obj) <- "AUC"
obj@meta.data$orig.ident=factor(obj@meta.data$orig.ident,levels =c(sample1,sample2))

# 作差异分析
library(limma)
data=obj@assays$AUC@data
dd_CM<- as.data.frame(as.matrix(data))
colnames(dd_CM)
table(Idents(obj))
group2<-  factor(c(rep('Je', 10107),rep('Se', 8838)))#顺序靠后的Se vs 顺序靠前的Je
dgelist <- DGEList(counts =dd_CM, group = group2)
#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
#（3）标准化，以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#差异表达基因分析
design <- model.matrix(~group2)
#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
write.table(lrt, sprintf('%svs%s_1.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)

# 根据差异分析结果选择差异明显的通路进行作图
features=c("BIOCARTA-VOBESITY-PATHWAY","WANG-CLASSIC-ADIPOGENIC-TARGETS-OF-PPARG","VERNOCHET-ADIPOGENESIS","GOBP-FAT-PAD-DEVELOPMENT","GOBP-NEGATIVE-REGULATION-OF-BROWN-FAT-CELL-DIFFERENTIATION")
Idents(obj)='orig.ident'
pdf(file=sprintf("%s-ActivepathwayFC0.4.pdf",obj),height = 6, width = 6)
DotPlot(obj,group.by = 'orig.ident', features = unique(features)) + coord_flip()+labs(x="Pathways")
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()


# 展示特定基因Nucb2在各类群中的表达情况
pdf(file="scRNASiJiAP.subanno-Nucb2-clusterdot.pdf",height = 4, width = 6)
DotPlot(scRNASiJiAP,group.by = 'seurat_clusters', features = unique(c('Nucb2'))) + coord_flip()+RotatedAxis()
dev.off()

# 不同基因之间的相关性展示
Idents(scRNASiJimyoFB) <- "cellType"
pdf(file="scRNASiJimyoFB_SiJi_Nucb2Igf1r.cor.pdf",height = 4, width = 6)
FeatureScatter(scRNASiJimyoFB, feature1 = "Nucb2", feature2 = "Igf1r")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
dev.off()
Idents(scRNASiJiEC) <- "cellType"
pdf(file="scRNASiJiEC_SiJi_Nucb2Vegfa.cor.pdf",height = 4, width = 6)
FeatureScatter(scRNASiJiEC, feature1 = "Nucb2", feature2 = "Vegfa")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
dev.off()





# 查看细胞数量
table(Idents(scRNASiJiAP))

# 不同基因表达的细胞类型占比
# EC总细胞数量
table(scRNASiJiEC@meta.data$cellType)
table(Idents(scRNASiJiEC))
filtered_cells <- WhichCells(scRNASiJiEC, expression = Nucb2 > 0)
length(filtered_cells)
filtered_cells <- WhichCells(scRNASiJiEC, expression = Vegfa > 0)
length(filtered_cells)
filtered_cells <- WhichCells(scRNASiJiEC, expression = Vegfa > 0&Nucb2 > 0)
length(filtered_cells)
filtered_cells <- WhichCells(scRNASiJiEC, expression = Vegfa == 0&Nucb2 == 0)
length(filtered_cells)
# 2433 
# 113 224 12 2108
# 101 212 12 2108
df=data.frame(name=c('EC_VEGF','EC_VEGF','EC_VEGF','EC_VEGF'),Type=c('Nucb2_s','Vegfa_s','DoubleT','DoubleF'),Num=c(101,212,12,2108))
data3 = ddply(df,'name',transform,percent_Num=Num/sum(Num)*100)
p <- ggplot(data3,aes(x=name,y=percent_Num,fill=Type))+geom_bar(stat = 'identity',width = 0.5,colour='black')+theme_classic()+labs(x='',y='Percentage')+theme(axis.title = element_text(size=20),axis.text = element_text(size=20))+scale_y_continuous(breaks=seq(0,100,25),labels=c('0','25%','50%','75%','100%'))+geom_text(aes(label=paste0(sprintf("%.1f", percent_Num), "%")),position = position_stack(vjust = 0.5), size = 8,color="black")+scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72"))+theme( legend.text = element_text(size = 18))
ggsave(filename = "scRNASiJiEC_VEGF_cor.ratio.pdf", plot = p, height = 6, width = 6)


# myoFB总细胞数量
table(scRNASiJimyoFB@meta.data$cellType)
table(Idents(scRNASiJimyoFB))
filtered_cells <- WhichCells(scRNASiJimyoFB, expression = Nucb2 > 0)
length(filtered_cells)
filtered_cells <- WhichCells(scRNASiJimyoFB, expression = Igf1r > 0)
length(filtered_cells)
filtered_cells <- WhichCells(scRNASiJimyoFB, expression = Igf1r > 0&Nucb2 > 0)
length(filtered_cells)
filtered_cells <- WhichCells(scRNASiJimyoFB, expression = Igf1r == 0&Nucb2 == 0)
length(filtered_cells)
# 17 
# 1 7 0 9 
# 1 7 0 9
df=data.frame(name=c('myoFB_IGF','myoFB_IGF','myoFB_IGF','myoFB_IGF'),Type=c('Nucb2_s','Igf1r_s','DoubleT','DoubleF'),Num=c(1,7,0,9))
data3 = ddply(df,'name',transform,percent_Num=Num/sum(Num)*100)
p <- ggplot(data3,aes(x=name,y=percent_Num,fill=Type))+geom_bar(stat = 'identity',width = 0.5,colour='black')+theme_classic()+labs(x='',y='Percentage')+theme(axis.title = element_text(size=20),axis.text = element_text(size=20))+scale_y_continuous(breaks=seq(0,100,25),labels=c('0','25%','50%','75%','100%'))+geom_text(aes(label=paste0(sprintf("%.1f", percent_Num), "%")),position = position_stack(vjust = 0.5), size = 8,color="black")+scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72"))+theme( legend.text = element_text(size = 18))
ggsave(filename = "scRNASiJimyoFB_IGF_cor.ratio.pdf", plot = p, height = 6, width = 6)








# ################4samples样本 ################
dir = c(
  "../Je/filtered_feature_bc_matrix",
  "../Ji/filtered_feature_bc_matrix",
  "../Se/filtered_feature_bc_matrix",
  "../Si/filtered_feature_bc_matrix"
)
scRNAlist <- list()  # 创建一个空的列表,其中包含三个seurat对象
names(dir) = c("Je", "Ji",'Se','Si')
sample_name = c("Je", "Ji",'Se','Si')
for (i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = sample_name[i])
}

# 合并Seurat object
scRNA <- merge(scRNAlist[[1]],y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]]))
dim(scRNA)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
# Visualize QC metrics as a violin plot
p_QC <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "p_QC.pdf", height = 6, width = 6, plot = p_QC)
saveRDS(file='scRNA_rawdata4.rds',scRNA)

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
saveRDS(scRNA,"scRNA4.filt.anchors.rds")

scRNA4=scRNA
scRNA=readRDS('../scRNA_seurat_cellType.RDS')
scRNASeJeAP=scRNA
SeuratObject=scRNA4

# scRNASeJeAP@meta.data$cellType2<- paste0(scRNASeJeAP@meta.data$cellType,scRNASeJeAP@meta.data$seurat_clusters)
scRNASeJeAP@meta.data$cellType2<-''
# 将df1转换为tibble，以便使用left_join，并设置行名为一个临时列
scRNASeJeAP_tibble <- tibble::rownames_to_column(as_tibble(scRNASeJeAP@meta.data), var = "rowname")
scRNASeJeAP_tibble$rowname=rownames(scRNASeJeAP@meta.data)
SeuratObject_tibble <- tibble::rownames_to_column(as_tibble(SeuratObject@meta.data), var = "rowname")
SeuratObject_tibble$rowname=rownames(SeuratObject@meta.data)
# 使用left_join合并数据框，然后根据结果替换df2中的celltype2列
merged_df <- left_join(SeuratObject_tibble, scRNASeJeAP_tibble, by = "rowname")
 
# 检查是否有NA值，如果有，则用原始的celltype2值填充
# merged_df$cellType <- as.character(merged_df$cellType)
# merged_df$cellType.x<- as.character(merged_df$cellType.x)
# # 使用ifelse进行替换
# merged_df$celltype2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)
 
# 如果需要，将celltype2转换回因子（并更新水平集）
# 注意：这一步只在celltype2原本就是因子且你需要保持它为因子类型时才需要
levels_to_use <- unique(c(as.character(merged_df$celltype2), levels(as.factor(scRNASeJeAP@meta.data$cellType))))
merged_df$cellType <- factor(merged_df$cellType, levels = levels_to_use)
SeuratObject@meta.data$cellType=merged_df$cellType
saveRDS(SeuratObject,"scRNA4.celltype.rds")

pdf(file="p_cellType_anno.pdf",height = 6, width = 6)
DimPlot(SeuratObject,group.by = "cellType",label=T)
dev.off()

pdf(file="p_seurat_clusters.pdf",height = 6, width = 6)
DimPlot(SeuratObject,group.by = "seurat_clusters",label=F)
dev.off()

DefaultAssay(SeuratObject) <- "RNA"
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

SeuratObject@meta.data$cellType=factor(SeuratObject@meta.data$cellType,levels =c("AP","B","mDC","EC","FB","MP","Mono","muscle","NK","myoFB","SMC","T","Unknown"))
pdf(file="p_genedotplot-cellType.pdf",height = 12, width = 9.75)
DotPlot(SeuratObject,group.by = 'cellType', features = unique(genes_to_check))+ coord_flip() +RotatedAxis()+ theme(text = element_text(size = 16))+theme(axis.text.x = element_text(size=20))+theme(axis.text.y = element_text(size=20))
dev.off()


####各细胞类型占比####
library(Seurat)
library(ggplot2)
library(dplyr)
SeuratObject=readRDS("scRNA4.celltype.rds")
scRNA=SeuratObject
# scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c("Se","Je","Si","Ji","ST","JT"))
# scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c("Se","Si","ST","Je","Ji","JT"))
scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c("Je","Se","Ji","Si"))
Ratio <- scRNA@meta.data %>%group_by(orig.ident,cellType) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(Freq = n/sum(n)*100)
write.csv(file='ratio1.csv',Ratio)

# 添加label
p1=ggplot(Ratio, aes(x = orig.ident, y = Freq, fill = cellType))+
  geom_col()+
  geom_text(aes(label = paste(round(Freq, 1),"%")), 
            position = position_stack(vjust = 0.5))+
  theme_classic()+
  scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                              "#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                              "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02"))
ggsave(filename = "p_Ratio.pdf", plot = p1, height = 6, width = 6)

# 不添加label
p1=ggplot(Ratio, aes(x = orig.ident, y = Freq, fill = cellType))+
  geom_col()+
  theme_classic()+
  scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                              "#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                              "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02"))
ggsave(filename = "p_Rationolabel.pdf", plot = p1, height = 6, width = 6)


#####按样本拆分scRNA#####
scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))
scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
# scRNASTJT=subset(x = scRNA, subset = (orig.ident == "ST"|orig.ident == "JT"))

pdf(file="Nucb2_scRNASeJe.pdf",height = 6, width = 12)
FeaturePlot(object = scRNASeJe, features = "Nucb2",split.by = 'orig.ident')
dev.off()

pdf(file="Nucb2_scRNASiJi.pdf",height = 6, width = 12)
FeaturePlot(object = scRNASiJi, features = "Nucb2",split.by = 'orig.ident')
dev.off()

# pdf(file="Nucb2_scRNASTJT.pdf",height = 6, width =12)
# FeaturePlot(object = scRNASTJT, features = "Nucb2",split.by = 'orig.ident')
# dev.off()





#合并两样本dimplot,展示gene
p1=DimPlot(scRNASeJe,group.by = 'orig.ident',
           cols = c("grey","grey"),
           pt.size = 1)+NoLegend()+ labs(title="Nucb2")

gene1 <- "Nucb2"

pos_gene=scRNASeJe@reductions$umap@cell.embeddings

pos_gene1=pos_gene[scRNASeJe@assays$RNA@data[gene1,]>0,]

pos_gene1Se=pos_gene1[grep('Se',rownames(pos_gene1)),]
pos_gene1Je=pos_gene1[grep('Je',rownames(pos_gene1)),]

p=p1+geom_point(data = as.data.frame(pos_gene1Se),
              aes(x=UMAP_1,y=UMAP_2), 
              shape = 16,size = 1.4, color="#6a5acd")+
  geom_point(data = as.data.frame(pos_gene1Je),
             aes(x=UMAP_1,y=UMAP_2),
             shape = 16,size = 1.4, color="#ff69b4")+theme(text = element_text(size = 30))+theme(axis.text.x = element_text(size=20),
           axis.text.y = element_text(size=20))
ggsave(filename = "SeJe-Dimplotgene.pdf", height = 6, width = 6, plot =p)

#合并两样本dimplot,展示gene
p1=DimPlot(scRNASiJi,group.by = 'orig.ident',
           cols = c("grey","grey"),
           pt.size = 1)+NoLegend()+ labs(title="Nucb2")

gene1 <- "Nucb2"

pos_gene=scRNASiJi@reductions$umap@cell.embeddings

pos_gene1=pos_gene[scRNASiJi@assays$RNA@data[gene1,]>0,]

pos_gene1Se=pos_gene1[grep('Si',rownames(pos_gene1)),]
pos_gene1Je=pos_gene1[grep('Ji',rownames(pos_gene1)),]

p=p1+geom_point(data = as.data.frame(pos_gene1Se),
              aes(x=UMAP_1,y=UMAP_2), 
              shape = 16,size = 1.4, color="#6a5acd")+
  geom_point(data = as.data.frame(pos_gene1Je),
             aes(x=UMAP_1,y=UMAP_2),
             shape = 16,size = 1.4, color="#ff69b4")+theme(text = element_text(size = 30))+theme(axis.text.x = element_text(size=20),
           axis.text.y = element_text(size=20))
ggsave(filename = "SiJi-Dimplotgene.pdf", height = 6, width = 6, plot =p)


# 通路活性分析
library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(SeuratData)

fatty <- read.gmt("fat.gmt")
head(fatty$term)
head(fatty$gene)
geneSets <- lapply(unique(fatty$term), function(x){print(x);fatty$gene[fatty$term == x]})
names(geneSets) <- unique(fatty$term)
features=rownames(unique(fatty$term))
features=gsub( "_", "-",features)#后续要加入到seurat，行名不能含有下划线

obj=scRNASeJe # 根据不同的对象进行调整
sample1='Je'
sample2='Se'

cells_rankings <- AUCell_buildRankings(obj@assays$RNA@data, plotStats=TRUE) 
# cells_rankings
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
aucs <- getAUC(cells_AUC)
# 将aucs加入seurat
# aucs=t(aucs)
obj_AUC <- CreateAssayObject(counts = aucs)
obj@assays$AUC <- obj_AUC
obj@assays$AUC@key <- "rna_"
DefaultAssay(obj) <- "AUC"
obj@meta.data$orig.ident=factor(obj@meta.data$orig.ident,levels =c(sample1,sample2))

# 作差异分析
library(limma)
data=obj@assays$AUC@data
dd_CM<- as.data.frame(as.matrix(data))
colnames(dd_CM)
table(Idents(obj))
group2<-  factor(c(rep('Je', 10107),rep('Se', 8838)))#顺序靠后的Se vs 顺序靠前的Je
dgelist <- DGEList(counts =dd_CM, group = group2)
#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
#（3）标准化，以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#差异表达基因分析
design <- model.matrix(~group2)
#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
write.table(lrt, sprintf('%svs%s_1.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)

# 根据差异分析结果选择差异明显的通路进行作图
features=c("BIOCARTA-VOBESITY-PATHWAY","WANG-CLASSIC-ADIPOGENIC-TARGETS-OF-PPARG","VERNOCHET-ADIPOGENESIS","GOBP-FAT-PAD-DEVELOPMENT","GOBP-NEGATIVE-REGULATION-OF-BROWN-FAT-CELL-DIFFERENTIATION")
Idents(obj)='orig.ident'
pdf(file=sprintf("%s-ActivepathwayFC0.4.pdf",obj),height = 6, width = 6)
DotPlot(obj,group.by = 'orig.ident', features = unique(features)) + coord_flip()+labs(x="Pathways")
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()


# 展示特定基因Nucb2在各类群中的表达情况
pdf(file="scRNASiJiAP.subanno-Nucb2-clusterdot.pdf",height = 4, width = 6)
DotPlot(scRNASiJiAP,group.by = 'seurat_clusters', features = unique(c('Nucb2'))) + coord_flip()+RotatedAxis()
dev.off()

# 不同基因之间的相关性展示
scRNASiJimyoFB=subset(x=scRNASiJi,subset=(cellType == "myoFB"))
scRNASiJiEC=subset(x=scRNASiJi,subset=(cellType == "EC"))

Idents(scRNASiJimyoFB) <- "cellType"
pdf(file="scRNASiJimyoFB_SiJi_Nucb2Igf1r.cor.pdf",height = 4, width = 6)
FeatureScatter(scRNASiJimyoFB, feature1 = "Nucb2", feature2 = "Igf1r")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
dev.off()
Idents(scRNASiJiEC) <- "cellType"
pdf(file="scRNASiJiEC_SiJi_Nucb2Vegfa.cor.pdf",height = 4, width = 6)
FeatureScatter(scRNASiJiEC, feature1 = "Nucb2", feature2 = "Vegfa")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
dev.off()





# 查看细胞数量
table(Idents(scRNASiJiAP))

# 不同基因表达的细胞类型占比
# EC总细胞数量
table(scRNASiJiEC@meta.data$cellType)
table(Idents(scRNASiJiEC))
filtered_cells <- WhichCells(scRNASiJiEC, expression = Nucb2 > 0)
length(filtered_cells)
filtered_cells <- WhichCells(scRNASiJiEC, expression = Vegfa > 0)
length(filtered_cells)
filtered_cells <- WhichCells(scRNASiJiEC, expression = Vegfa > 0&Nucb2 > 0)
length(filtered_cells)
filtered_cells <- WhichCells(scRNASiJiEC, expression = Vegfa == 0&Nucb2 == 0)
length(filtered_cells)
# 2433 
# 113 224 12 2108
# 101 212 12 2108
df=data.frame(name=c('EC_VEGF','EC_VEGF','EC_VEGF','EC_VEGF'),Type=c('Nucb2_s','Vegfa_s','DoubleT','DoubleF'),Num=c(101,212,12,2108))
data3 = ddply(df,'name',transform,percent_Num=Num/sum(Num)*100)
p <- ggplot(data3,aes(x=name,y=percent_Num,fill=Type))+geom_bar(stat = 'identity',width = 0.5,colour='black')+theme_classic()+labs(x='',y='Percentage')+theme(axis.title = element_text(size=20),axis.text = element_text(size=20))+scale_y_continuous(breaks=seq(0,100,25),labels=c('0','25%','50%','75%','100%'))+geom_text(aes(label=paste0(sprintf("%.1f", percent_Num), "%")),position = position_stack(vjust = 0.5), size = 8,color="black")+scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72"))+theme( legend.text = element_text(size = 18))
ggsave(filename = "scRNASiJiEC_VEGF_cor.ratio.pdf", plot = p, height = 6, width = 6)


# myoFB总细胞数量
table(scRNASiJimyoFB@meta.data$cellType)
table(Idents(scRNASiJimyoFB))
filtered_cells <- WhichCells(scRNASiJimyoFB, expression = Nucb2 > 0)
length(filtered_cells)
filtered_cells <- WhichCells(scRNASiJimyoFB, expression = Igf1r > 0)
length(filtered_cells)
filtered_cells <- WhichCells(scRNASiJimyoFB, expression = Igf1r > 0&Nucb2 > 0)
length(filtered_cells)
filtered_cells <- WhichCells(scRNASiJimyoFB, expression = Igf1r == 0&Nucb2 == 0)
length(filtered_cells)
# 17 
# 1 7 0 9 
# 1 7 0 9
df=data.frame(name=c('myoFB_IGF','myoFB_IGF','myoFB_IGF','myoFB_IGF'),Type=c('Nucb2_s','Igf1r_s','DoubleT','DoubleF'),Num=c(1,7,0,9))
data3 = ddply(df,'name',transform,percent_Num=Num/sum(Num)*100)
p <- ggplot(data3,aes(x=name,y=percent_Num,fill=Type))+geom_bar(stat = 'identity',width = 0.5,colour='black')+theme_classic()+labs(x='',y='Percentage')+theme(axis.title = element_text(size=20),axis.text = element_text(size=20))+scale_y_continuous(breaks=seq(0,100,25),labels=c('0','25%','50%','75%','100%'))+geom_text(aes(label=paste0(sprintf("%.1f", percent_Num), "%")),position = position_stack(vjust = 0.5), size = 8,color="black")+scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72"))+theme( legend.text = element_text(size = 18))
ggsave(filename = "scRNASiJimyoFB_IGF_cor.ratio.pdf", plot = p, height = 6, width = 6)


# 相互作用分析
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)
searchsabbr=c('SS','ECM-R','CCC')
# 比较分组的样本名称，用于保存文件
comparegroup=c('SeJe','SiJi','STJT')
# 比较分组的样本在object.list所在的位置进行分组，用于后续分析
groups=list(c(1,4),c(2,5),c(3,6))
# searchs=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
searchs=c('Secreted Signaling')
for (i in 1:length(searchs)){
    search=searchs[i]
    abbr=searchsabbr[i]
    # search=searchs[1]
    # abbr=searchsabbr[1]
    # 整合Se Si ST Je Ji JT Secreted Signaling cellchat
    load(sprintf('../cellchat/Se_cellchat_%s.RData',search))
    Se_cellchat=cellchat
    load(sprintf('../cellchat/Si_cellchat_%s.RData',search))
    Si_cellchat=cellchat
    load(sprintf('../cellchat/ST_cellchat_%s.RData',search))
    ST_cellchat=cellchat
    load(sprintf('../cellchat/Je_cellchat_%s.RData',search))
    Je_cellchat=cellchat
    load(sprintf('../cellchat/Ji_cellchat_%s.RData',search))
    Ji_cellchat=cellchat
    load(sprintf('../cellchat/JT_cellchat_%s.RData',search))
    JT_cellchat=cellchat

    object.list <- list(Je=Je_cellchat,Ji=Ji_cellchat,JT=JT_cellchat,Se=Se_cellchat,Si =Si_cellchat,ST=ST_cellchat)
    #run netAnalysis_computeCentrality
    object.list<- lapply(object.list,function(x){
        x=netAnalysis_computeCentrality(x)})
    cellchat_m <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
    ####从宏观角度预测细胞通讯####
    #比较交互总数和交互强度

    for (num in 1:length(groups)){
        compare=comparegroup[num]
        group=groups[[num]]
        # compare=comparegroup[3]
        # group=groups[[3]]
        pdf(file=sprintf("%sChatNumall_%s.pdf",compare,abbr),height = 8, width = 8)
        par(mfrow = c(1,1), xpd=TRUE)
        compareInteractions(cellchat_m, show.legend = F, group = (1:6),size.text = 20)#group颜色向量 默认measure='count'
        print(compareInteractions(cellchat_m, show.legend = F, group = (1:6),size.text = 20))
        compareInteractions(cellchat_m, show.legend = F, group = (1:6), measure = "weight",size.text = 20)
        print(compareInteractions(cellchat_m, show.legend = F, group = (1:6), measure = "weight",size.text = 20))
        dev.off()
        #不同细胞群之间的相互作用数量或强度的差异 circle
        pdf(file=sprintf("%sChatNumCell_%s.pdf",compare,abbr),height = 8, width = 8)
        par(mfrow = c(1,1), xpd=TRUE)
        netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group)#group为相互比较的不同的cellchat，通过names(object.list)查看用来比较的cellchat所在的位置
        # print(netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group))
        netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group)
        # print(netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group))
        dev.off() 
        #不同细胞群之间的相互作用数量或强度的差异 circle
        pdf(file=sprintf("%sChatNumCellDIFF-AP_%s.pdf",compare,abbr),height = 6, width = 6)
        par(mfrow = c(1,1), xpd=TRUE)
        netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group,sources.use='AP')#group为相互比较的不同的cellchat，通过names(object.list)查看用来比较的cellchat所在的位置,compare[2]-compare[1]
        # print(netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group))
        netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group,sources.use='AP')
        # print(netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group))

        netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group,targets.use='AP')#group为相互比较的不同的cellchat，通过names(object.list)查看用来比较的cellchat所在的位置
        # print(netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group))
        netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group,targets.use='AP')
        # print(netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group))
        dev.off()
        #bubble        
        pdf(file=sprintf("%sChatbubble_%s-AME.pdf",compare,abbr),height = 8, width = 8)
        par(mfrow = c(1,1), xpd=TRUE)
        netVisual_bubble(cellchat_m, sources.use =c('AP','myoFB','EC'),targets.use = c('AP','myoFB','EC'),comparison = group, angle.x = 45, font.size = 16) +
        theme(legend.title = element_text(size = 16), legend.text = element_text(size = 16))
        print(netVisual_bubble(cellchat_m, sources.use =c('AP','myoFB','EC'),targets.use = c('AP','myoFB','EC'),comparison =  group, angle.x = 45, font.size = 16)) +
        theme(legend.title = element_text(size = 16), legend.text = element_text(size = 16))
        dev.off()
      }
}



# gsva分析
library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
library(patchwork)
library(clusterProfiler)
library(GSVA)
scRNA <- readRDS("../scRNA_seurat_cellType.RDS")
scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
genelist<- list.files('../cellchat/gmtps/')

#Mm_m5 <- msigdbr(species = "Mus musculus", category = "m5")## 定义基因集，选取C2
#fgsea_sets = mdb_c2 [grep("^GO",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets<- list()
for (i in 1:length(genelist)) {
  dir<- paste0('../cellchat/gmtps/',genelist[i])
  temp<- read.gmt(dir)
  fgsea_sets[[names(table(temp$term))]]<- temp$gene
  
}
scRNASiJimyoFB=subset(x=scRNASiJi,subset=(cellType == "myoFB"))
data=scRNASiJimyoFB@assays$RNA@data

expr=as.matrix(data) 
kegg <- gsva(expr, fgsea_sets, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
write.csv(file='kegg.csv',kegg)
# p=pheatmap(kegg)#绘制热图
Idents(scRNASiJimyoFB)='orig.ident'
CellsClusters <- data.frame(Cell = names(Idents(scRNASiJimyoFB)), 
                            CellType = as.character(Idents(scRNASiJimyoFB)),
                            stringsAsFactors = FALSE)
kegg_scores_df <- as.data.frame(t(kegg)) %>% rownames_to_column("Cell") %>% gather(Pathway, Activity, -Cell)
kegg_scores_df <- inner_join(kegg_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- kegg_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
#########plot the different pathway activities for the different cell populations
## We prepare the data for the plot
summarized_kegg_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_kegg_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_kegg_scores_df)/paletteLength, 
                      max(summarized_kegg_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_kegg_scores_df),fontsize=16, 
                        fontsize_row = 16, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "GSVA", angle_col = 0,
                        treeheight_col = 0,  border_color = NA)
ggsave('SiJi-myoFB-gsva_hmap.pdf',progeny_hmap,height=6,width=10)

progeny_hmap = pheatmap(t(summarized_kegg_scores_df),fontsize=16, 
                        fontsize_row = 16, 
                        color=myColor, 
                        main = "GSVA", angle_col = 0,
                        treeheight_col = 0,  border_color = NA)
ggsave('SiJi-myoFB-gsva_hmap2.pdf',progeny_hmap,height=6,width=10)

# GSVA
genelist<- list.files('../cellchat/gmtps/')

#Mm_m5 <- msigdbr(species = "Mus musculus", category = "m5")## 定义基因集，选取C2
#fgsea_sets = mdb_c2 [grep("^GO",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets<- list()
for (i in 1:length(genelist)) {
  dir<- paste0('../cellchat/gmtps/',genelist[i])
  temp<- read.gmt(dir)
  fgsea_sets[[names(table(temp$term))]]<- temp$gene
  
}
scRNASiJiEC=subset(x=scRNASiJi,subset=(cellType == "EC"))
data=scRNASiJiEC@assays$RNA@data

expr=as.matrix(data) 
kegg <- gsva(expr, fgsea_sets, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
write.csv(file='kegg.csv',kegg)
# p=pheatmap(kegg)#绘制热图
Idents(scRNASiJiEC)='orig.ident'
CellsClusters <- data.frame(Cell = names(Idents(scRNASiJiEC)), 
                            CellType = as.character(Idents(scRNASiJiEC)),
                            stringsAsFactors = FALSE)
kegg_scores_df <- as.data.frame(t(kegg)) %>% rownames_to_column("Cell") %>% gather(Pathway, Activity, -Cell)
kegg_scores_df <- inner_join(kegg_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- kegg_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
#########plot the different pathway activities for the different cell populations
## We prepare the data for the plot
summarized_kegg_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_kegg_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_kegg_scores_df)/paletteLength, 
                      max(summarized_kegg_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_kegg_scores_df),fontsize=16, 
                        fontsize_row = 16, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "GSVA", angle_col = 0,
                        treeheight_col = 0,  border_color = NA)
ggsave('SiJi-EC-gsva_hmap.pdf',progeny_hmap,height=6,width=10)

progeny_hmap = pheatmap(t(summarized_kegg_scores_df),fontsize=16, 
                        fontsize_row = 16, 
                        color=myColor,  
                        main = "GSVA", angle_col = 0,
                        treeheight_col = 0,  border_color = NA)
ggsave('SiJi-EC-gsva_hmap2.pdf',progeny_hmap,height=6,width=10)

# 差异基因分析
SeJemyoFB Nucb2
SeJeAP
SiJimyoFB Nucb2
SiJiAP

scRNA <- readRDS("../scRNA_seurat_cellType.RDS")

scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))

SeJemyoFB=subset(x = scRNASeJe, subset = (cellType == "myoFB"))
obj=SeJemyoFB
Idents(obj)='orig.ident'
table(Idents(obj))
Se  4817  Je  3

SeJeAP=subset(x = scRNASeJe, subset = (cellType == "AP"))
obj=SeJeAP
Idents(obj)='orig.ident'
table(Idents(obj))
 Se   Je 
 612 4162

scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
SiJimyoFB=subset(x = scRNASiJi, subset = (cellType == "myoFB"))
obj=SiJimyoFB
Idents(obj)='orig.ident'
table(Idents(obj))
Si Ji 
16  1 


SiJiAP=subset(x = scRNASiJi, subset = (cellType == "AP"))
obj=SiJiAP
Idents(obj)='orig.ident'
table(Idents(obj))
  Si   Ji 
 526 5069 


library(limma)
data=obj@assays$RNA@data
dd_CM<- as.data.frame(as.matrix(data))
# colnames(dd_CM)

group2<-  factor(c(rep('Je', 10107),rep('Se', 8838)))#顺序靠后的Se vs 顺序靠前的Je
dgelist <- DGEList(counts =dd_CM, group = group2)
#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
#（3）标准化，以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#差异表达基因分析
design <- model.matrix(~group2)
#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
write.table(lrt, sprintf('%svs%s_1.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)
library(ggplot2)

d1<- read.table('RHvsM_2.txt')
colnames(d1)
d1$gene<- rownames(d1)

# 设置p_value和logFC的阈值
cut_off_FDR = 0.05  #统计显著性
cut_off_logFC = 0.5           #差异倍数值
d1$change<- 'Stable'
# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
d1$change = ifelse(d1$FDR< cut_off_FDR & abs(d1$logFC) > cut_off_logFC, 
                        ifelse(d1$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(d1)
table(d1$change)

p <- ggplot(
  # 数据、映射、颜色
  d1, aes(x = logFC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_FDR),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  scale_x_continuous(limits = c(-5, 5))+
  scale_y_continuous(limits = c(0,30))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('RH v.s. M')
# p

# 将需要标记的基因放置在label列(logFC >= 5)
library(ggrepel)
d1$label <- ifelse(d1$PValue < cut_off_FDR & abs(d1$logFC) >= 0.3,
                        as.character(d1$gene), "")


p=p + geom_label_repel(data = d1, aes(x = d1$logFC, 
                                         y = -log10(d1$PValue), 
                                         label = label),
                     size = 3, box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"), 
                     segment.color = "black", 
                     show.legend = FALSE)
ggsave(filename = "diffgene.pdf", height = 6, width = 6, plot = p)





library(monocle)
library(RColorBrewer)
library(Seurat)
library(tidyverse)
library(patchwork)
# AP myoFB细胞分化轨迹
scRNA <- readRDS("../scRNA_seurat_cellType.RDS")
scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))
scRNASeJeAPmyoFB=subset(x = scRNASeJe, subset = (cellType == "AP"|cellType == "myoFB"))
SeuratObject=scRNASeJeAPmyoFB
scRNASeJeAP=readRDS('../cellchat/scRNASeJeAPharmony.rda')

# SeuratObject$Cluster<- paste0('cluster',SeuratObject$seurat_clusters)
scRNASeJeAP@meta.data$cellType2<- paste0(scRNASeJeAP@meta.data$cellType,scRNASeJeAP@meta.data$seurat_clusters)
library(dplyr)
 
# 将df1转换为tibble，以便使用left_join，并设置行名为一个临时列
scRNASeJeAP_tibble <- tibble::rownames_to_column(as_tibble(scRNASeJeAP@meta.data), var = "rowname")
scRNASeJeAP_tibble$rowname=rownames(scRNASeJeAP@meta.data)
SeuratObject_tibble <- tibble::rownames_to_column(as_tibble(SeuratObject@meta.data), var = "rowname")
SeuratObject_tibble$rowname=rownames(SeuratObject@meta.data)
# 使用left_join合并数据框，然后根据结果替换df2中的celltype2列
merged_df <- left_join(SeuratObject_tibble, scRNASeJeAP_tibble, by = "rowname")
 
# 检查是否有NA值，如果有，则用原始的celltype2值填充
merged_df$cellType2 <- as.character(merged_df$cellType2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
# 使用ifelse进行替换
merged_df$celltype2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)
 
# 如果需要，将celltype2转换回因子（并更新水平集）
# 注意：这一步只在celltype2原本就是因子且你需要保持它为因子类型时才需要
levels_to_use <- unique(c(as.character(merged_df$celltype2), levels(as.factor(scRNASeJeAP@meta.data$cellType))))
merged_df$celltype2 <- factor(merged_df$celltype2, levels = levels_to_use)
SeuratObject@meta.data$celltype2=merged_df$celltype2

expr_matrix <- as.data.frame(SeuratObject@assays$RNA@counts)
# 细胞注释矩阵（列为细胞名）
sample_sheet <- SeuratObject@meta.data
# 基因注释矩阵（行为基因名）
gene_annotation <- data.frame(
  gene_short_name=rownames(SeuratObject@assays$RNA),
  row.names = rownames(SeuratObject@assays$RNA)
)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

#创建monocle对象
cd <- newCellDataSet(as(as.matrix(expr_matrix),"sparseMatrix"), 
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

###计算数据的size factors和dispersions
cd <- estimateSizeFactors(cd)
cd <- estimateDispersions(cd)
# 过滤低于1%细胞中检出的基因，最低表达阈值为0.5
cd <- detectGenes(cd, min_expr = 0.5)
expressed_genes <- row.names(subset(fData(cd), num_cells_expressed > nrow(sample_sheet) * 0.01))
# 查看筛选后的基因个数（用于后面基因的筛选）
length(expressed_genes)

###构建轨迹
disp_table <- dispersionTable(cd)
# 高度离散基因的筛选标准，可根据数据情况设置mean_expression的值
ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)
# 根据数据的某种分类进行差异分析
diff_test_res <- differentialGeneTest(cd[expressed_genes,],fullModelFormulaStr = "~celltype2")
# 筛选差异基因（q < 1e-5并且属于之前计算的expressed_genes列表中）
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-5))
ordering_genes <- intersect(ordering_genes, expressed_genes)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)
# 选取的基因数目为每个细胞的维度，基于默认的'DDRTree'方法进行数据降维
cd <- reduceDimension(cd, max_components = 2, method = 'DDRTree')
# 对细胞进行排序，由于排序无法区分起点和终点，若分析所得时序与实际相反，根据“reverse”参数进行调整，默认reverse=F
cd <- orderCells(cd, reverse = F) 
# 排序好的细胞可以进行可视化，可标注细胞的各注释信息)
# 查看细胞注释信息
head(cd@phenoData@data)
saveRDS(cd,file='SeJem2_cd_apmyofb.rds')
cd <- readRDS('SeJem2_cd_apmyofb.rds')
# 设置颜色
color1 <- c(brewer.pal(8, "Set1"))
getPalette <- colorRampPalette(brewer.pal(6, "Set1"))
# # 不拆分
# # 按照表型进行映射
# pdf(file="cd_cellType2.trajectory.pdf",width=6.5,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "celltype2") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"celltype2"]))))
# dev.off()
# # 按照表型进行映射
# pdf(file="cd_cellType.trajectory.pdf",width=6.5,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"celltype2"]))))
# dev.off()

# 拆分
pData(cd)$sample=SeuratObject$orig.ident
pData(cd)$sample=factor(pData(cd)$sample, levels = c('Je','Se'))
pData(cd)$celltype2 <- factor(pData(cd)$celltype2, levels = c('AP4','AP2','AP1','AP7','AP0','AP5','AP3','AP6','myoFB'))

nrows_state=1
pdf(file="sejecd_cellType.trajectory.pdf",width=12,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"cellType"]))))
plot_cell_trajectory(cd, color_by = "cellType", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~sample, nrow = nrows_state, scales = "free")+ scale_color_manual(values=c('red','blue'))+ggtitle("cellType") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()

pdf(file="sejecd_sample.trajectory.pdf",width=12,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"cellType"]))))
plot_cell_trajectory(cd, color_by = "sample", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~celltype2, nrow = 3, scales = "free")+ggtitle("sample") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()



scRNASeJeAP=subset(x = scRNASeJe, subset = (cellType == "AP"))
SeuratObject=scRNASeJeAP
scRNASeJeAP=readRDS('../cellchat/scRNASeJeAPharmony.rda')

# SeuratObject$Cluster<- paste0('cluster',SeuratObject$seurat_clusters)
scRNASeJeAP@meta.data$cellType2<- paste0(scRNASeJeAP@meta.data$cellType,scRNASeJeAP@meta.data$seurat_clusters)
library(dplyr)
 
# 将df1转换为tibble，以便使用left_join，并设置行名为一个临时列
scRNASeJeAP_tibble <- tibble::rownames_to_column(as_tibble(scRNASeJeAP@meta.data), var = "rowname")
scRNASeJeAP_tibble$rowname=rownames(scRNASeJeAP@meta.data)
SeuratObject_tibble <- tibble::rownames_to_column(as_tibble(SeuratObject@meta.data), var = "rowname")
SeuratObject_tibble$rowname=rownames(SeuratObject@meta.data)
# 使用left_join合并数据框，然后根据结果替换df2中的celltype2列
merged_df <- left_join(SeuratObject_tibble, scRNASeJeAP_tibble, by = "rowname")
 
# 检查是否有NA值，如果有，则用原始的celltype2值填充
merged_df$cellType2 <- as.character(merged_df$cellType2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
# 使用ifelse进行替换
merged_df$celltype2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)
 
# 如果需要，将celltype2转换回因子（并更新水平集）
# 注意：这一步只在celltype2原本就是因子且你需要保持它为因子类型时才需要
levels_to_use <- unique(c(as.character(merged_df$celltype2), levels(as.factor(scRNASeJeAP@meta.data$cellType))))
merged_df$celltype2 <- factor(merged_df$celltype2, levels = levels_to_use)
SeuratObject@meta.data$celltype2=merged_df$celltype2

expr_matrix <- as.data.frame(SeuratObject@assays$RNA@counts)
# 细胞注释矩阵（列为细胞名）
sample_sheet <- SeuratObject@meta.data
# 基因注释矩阵（行为基因名）
gene_annotation <- data.frame(
  gene_short_name=rownames(SeuratObject@assays$RNA),
  row.names = rownames(SeuratObject@assays$RNA)
)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

#创建monocle对象
cd <- newCellDataSet(as(as.matrix(expr_matrix),"sparseMatrix"), 
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

###计算数据的size factors和dispersions
cd <- estimateSizeFactors(cd)
cd <- estimateDispersions(cd)
# 过滤低于1%细胞中检出的基因，最低表达阈值为0.5
cd <- detectGenes(cd, min_expr = 0.5)
expressed_genes <- row.names(subset(fData(cd), num_cells_expressed > nrow(sample_sheet) * 0.01))
# 查看筛选后的基因个数（用于后面基因的筛选）
length(expressed_genes)

###构建轨迹
disp_table <- dispersionTable(cd)
# 高度离散基因的筛选标准，可根据数据情况设置mean_expression的值
ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)
# 根据数据的某种分类进行差异分析
diff_test_res <- differentialGeneTest(cd[expressed_genes,],fullModelFormulaStr = "~celltype2")
# 筛选差异基因（q < 1e-5并且属于之前计算的expressed_genes列表中）
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-5))
ordering_genes <- intersect(ordering_genes, expressed_genes)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)
# 选取的基因数目为每个细胞的维度，基于默认的'DDRTree'方法进行数据降维
cd <- reduceDimension(cd, max_components = 2, method = 'DDRTree')
# 对细胞进行排序，由于排序无法区分起点和终点，若分析所得时序与实际相反，根据“reverse”参数进行调整，默认reverse=F
cd <- orderCells(cd, reverse = F) 
# 排序好的细胞可以进行可视化，可标注细胞的各注释信息)
# 查看细胞注释信息
head(cd@phenoData@data)
saveRDS(cd,file='sejem2_cd_ap.rds')
# 图如下
BEAM_res <- BEAM(cd, branch_point = 3, cores = 50, progenitor_method = "duplicate")
write.csv(file="SeJe.BEAM_res.csv",BEAM_res,rowname)
# BEAM_res=read.csv(file="SeJe.BEAM_res.csv")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
pdf(file="sejem2_genes_branched_heatmap50.pdf",width=4,height=6)
plot_genes_branched_heatmap(cd[c(row.names(subset(BEAM_res[1:50,],qval < 1e-4))),],
                            branch_point = 3,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()






# SiJiAPmyoFB
scRNA=readRDS('../scRNA_seurat_cellType.RDS')
scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
scRNASiJiAPmyoFB=subset(x=scRNASiJi,subset=(cellType == "AP"|cellType == "myoFB"))

SeuratObject=scRNASiJiAPmyoFB
scRNASiJiAP=readRDS('../APSIJI/scRNASiJiAPharmony.rda')


# object=scRNASiJiAP
# object$seurat_clusters <- as.character(object$seurat_clusters)
# object$seurat_clusters[object$seurat_clusters == "5"] <- "1"
# object$seurat_clusters[object$seurat_clusters == "8"] <- "5"
# object$seurat_clusters <- factor(object$seurat_clusters)



# SeuratObject$Cluster<- paste0('cluster',SeuratObject$seurat_clusters)
scRNASiJiAP@meta.data$cellType2<- paste0(scRNASiJiAP@meta.data$cellType,scRNASiJiAP@meta.data$seurat_clusters)
library(dplyr)
 
# 将df1转换为tibble，以便使用left_join，并设置行名为一个临时列
scRNASiJiAP_tibble <- tibble::rownames_to_column(as_tibble(scRNASiJiAP@meta.data), var = "rowname")
scRNASiJiAP_tibble$rowname=rownames(scRNASiJiAP@meta.data)
SeuratObject_tibble <- tibble::rownames_to_column(as_tibble(SeuratObject@meta.data), var = "rowname")
SeuratObject_tibble$rowname=rownames(SeuratObject@meta.data)
# 使用left_join合并数据框，然后根据结果替换df2中的celltype2列
merged_df <- left_join(SeuratObject_tibble, scRNASiJiAP_tibble, by = "rowname")
 
# 检查是否有NA值，如果有，则用原始的celltype2值填充
merged_df$cellType2 <- as.character(merged_df$cellType2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
# 使用ifelse进行替换
merged_df$celltype2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)
 
# 如果需要，将celltype2转换回因子（并更新水平集）
# 注意：这一步只在celltype2原本就是因子且你需要保持它为因子类型时才需要
levels_to_use <- unique(c(as.character(merged_df$celltype2), levels(as.factor(scRNASiJiAP@meta.data$cellType))))
merged_df$celltype2 <- factor(merged_df$celltype2, levels = levels_to_use)
SeuratObject@meta.data$celltype2=merged_df$celltype2

SeuratObject2=SeuratObject
DefaultAssay(SeuratObject2) <- "RNA"

SeuratObject2 <- NormalizeData(SeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

SeuratObject2 <- FindVariableFeatures(SeuratObject2,nfeatures = 3000)
SeuratObject2 <- ScaleData(SeuratObject2,vars.to.regress = 'percent.mt')
SeuratObject2 <- RunPCA(SeuratObject2)
SeuratObject2 <- FindNeighbors(SeuratObject2, reduction = "pca", dims = 1:30)
SeuratObject2 <- FindClusters(SeuratObject2, resolution = 0.3)
SeuratObject2 <- RunUMAP(SeuratObject2, reduction = "pca", dims = 1:30) 
library(harmony)
SeuratObject2 <- RunHarmony(SeuratObject2, group.by.vars = "orig.ident")
SeuratObject2  <- FindNeighbors(SeuratObject2 , dims = 1:30 , reduction = "harmony")
SeuratObject2  <- FindClusters(SeuratObject2, save.snn=T , resolution = 0.3)
SeuratObject2  <- RunUMAP(SeuratObject2 , dims=1:30,reduction='harmony')

p1=DimPlot(SeuratObject2,group.by = "orig.ident")
# 查看聚类信息
p2=DimPlot(SeuratObject2,group.by = "seurat_clusters",label=T)
p3=DimPlot(SeuratObject2,group.by = "celltype2",label=T)
combine<-CombinePlots(list(p1,p2,p3))
ggsave(filename = "sijip_DimPlot_sample_seurat_clusters.pdf", plot = combine, height = 6, width = 10)

library(monocle3)
SeuratObject=SeuratObject2
data <- GetAssayData(SeuratObject, assay = 'RNA', slot = 'counts')
cell_metadata <- SeuratObject@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 19)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP",color_cells_by="celltype2") + ggtitle('cds.umap')
########可视化部分基因
# ciliated_genes <- c("MITF",
#                     "SOX10",
#                     "AXL",
#                     "NGFR")

# plot_cells(cds,
#            genes=ciliated_genes,
#            label_cell_groups=FALSE,
#            show_trajectory_graph=FALSE)
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(SeuratObject, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype2") + ggtitle('int.umap')
ggsave(filename='sijicdsintumap.pdf',plot = p2, height = 6, width = 5)
# p1+p2
combine<-CombinePlots(list(p1,p2))
ggsave(filename = "sijicds_int.umap.pdf", plot = combine, height = 6, width = 10)

library(monocle)
detach("package:monocle3", unload = TRUE)

scRNA <- readRDS("../scRNA_seurat_cellType.RDS")
scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
scRNASiJiAPmyoFB=subset(x = scRNASiJi, subset = (cellType == "AP"|cellType == "myoFB"))
SeuratObject=scRNASiJiAPmyoFB
scRNASiJiAP=readRDS('')

# SeuratObject$Cluster<- paste0('cluster',SeuratObject$seurat_clusters)
scRNASiJiAP@meta.data$cellType2<- paste0(scRNASiJiAP@meta.data$cellType,scRNASiJiAP@meta.data$seurat_clusters)
library(dplyr)
 
# 将df1转换为tibble，以便使用left_join，并设置行名为一个临时列
scRNASiJiAP_tibble <- tibble::rownames_to_column(as_tibble(scRNASiJiAP@meta.data), var = "rowname")
scRNASiJiAP_tibble$rowname=rownames(scRNASiJiAP@meta.data)
SeuratObject_tibble <- tibble::rownames_to_column(as_tibble(SeuratObject@meta.data), var = "rowname")
SeuratObject_tibble$rowname=rownames(SeuratObject@meta.data)
# 使用left_join合并数据框，然后根据结果替换df2中的celltype2列
merged_df <- left_join(SeuratObject_tibble, scRNASiJiAP_tibble, by = "rowname")
 
# 检查是否有NA值，如果有，则用原始的celltype2值填充
merged_df$cellType2 <- as.character(merged_df$cellType2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
# 使用ifelse进行替换
merged_df$celltype2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)
 
# 如果需要，将celltype2转换回因子（并更新水平集）
# 注意：这一步只在celltype2原本就是因子且你需要保持它为因子类型时才需要
merged_df$celltype2 <- factor(merged_df$celltype2, levels = c('AP0','AP1','AP2','AP3','AP4','AP5','AP6','AP7','myoFB'))
SeuratObject@meta.data$celltype2=merged_df$celltype2
expr_matrix <- as.data.frame(SeuratObject@assays$RNA@counts)
# 细胞注释矩阵（列为细胞名）
sample_sheet <- SeuratObject@meta.data
# 基因注释矩阵（行为基因名）
gene_annotation <- data.frame(
  gene_short_name=rownames(SeuratObject@assays$RNA),
  row.names = rownames(SeuratObject@assays$RNA)
)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

#创建monocle对象
cd <- newCellDataSet(as(as.matrix(expr_matrix),"sparseMatrix"), 
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

###计算数据的size factors和dispersions
cd <- estimateSizeFactors(cd)
cd <- estimateDispersions(cd)
# 过滤低于1%细胞中检出的基因，最低表达阈值为0.5
cd <- detectGenes(cd, min_expr = 0.5)
expressed_genes <- row.names(subset(fData(cd), num_cells_expressed > nrow(sample_sheet) * 0.01))
# 查看筛选后的基因个数（用于后面基因的筛选）
length(expressed_genes)

###构建轨迹
disp_table <- dispersionTable(cd)
# 高度离散基因的筛选标准，可根据数据情况设置mean_expression的值
ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)
# 根据数据的某种分类进行差异分析
diff_test_res <- differentialGeneTest(cd[expressed_genes,],fullModelFormulaStr = "~celltype2")
# 筛选差异基因（q < 1e-5并且属于之前计算的expressed_genes列表中）
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-5))
ordering_genes <- intersect(ordering_genes, expressed_genes)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)
# 选取的基因数目为每个细胞的维度，基于默认的'DDRTree'方法进行数据降维
cd <- reduceDimension(cd, max_components = 2, method = 'DDRTree')
# 对细胞进行排序，由于排序无法区分起点和终点，若分析所得时序与实际相反，根据“reverse”参数进行调整，默认reverse=F
cd <- orderCells(cd, reverse = F) 
# 排序好的细胞可以进行可视化，可标注细胞的各注释信息)
# 查看细胞注释信息
head(cd@phenoData@data)
saveRDS(cd,file='SiJim2_cd_apmyofb.rds')
cd <- readRDS('SiJim2_cd_apmyofb.rds')
# 设置颜色
color1 <- c(brewer.pal(8, "Set1"))
getPalette <- colorRampPalette(brewer.pal(6, "Set1"))
# # 不拆分
# # 按照表型进行映射
# pdf(file="cd_cellType2.trajectory.pdf",width=6.5,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "celltype2") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"celltype2"]))))
# dev.off()
# # 按照表型进行映射
# pdf(file="cd_cellType.trajectory.pdf",width=6.5,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"celltype2"]))))
# dev.off()

cd=readRDS('sijim2_cd_apmyofb.rds')
# 拆分
pData(cd)$sample=SeuratObject$orig.ident
pData(cd)$sample=factor(pData(cd)$sample, levels = c('Ji','Si'))

pData(cd)$celltype2 <- as.character(pData(cd)$celltype2)
pData(cd)$celltype2[pData(cd)$celltype2 == "AP5"] <- "AP1"
pData(cd)$celltype2[pData(cd)$celltype2 == "AP8"] <- "AP5"
pData(cd)$celltype2 <- factor(pData(cd)$celltype2)
nrows_state=1
pdf(file="SiJicd_cellType.trajectory.pdf",width=12,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"cellType"]))))
plot_cell_trajectory(cd, color_by = "cellType", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~sample, nrow = nrows_state, scales = "free")+ scale_color_manual(values=c('red','blue'))+ggtitle("cellType") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()

0,2,6,8,1,3,4,5,7,myoFB
pData(cd)$celltype2 <- factor(pData(cd)$celltype2, levels = c('AP0','AP2','AP6','AP5','AP3','AP4','AP7','AP1','myoFB'))
pdf(file="SiJicd_sample.trajectory.pdf",width=12,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"cellType"]))))
plot_cell_trajectory(cd, color_by = "sample", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~celltype2, nrow = 3, scales = "free")+ggtitle("sample") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()




scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))
scRNASeJeAPmyoFB=subset(x = scRNASeJe, subset = (cellType == "AP"|cellType == "myoFB"))
SeuratObject=scRNASeJeAPmyoFB
scRNASeJeAP=readRDS('../cellchat/scRNASeJeAPharmony.rda')

# SeuratObject$Cluster<- paste0('cluster',SeuratObject$seurat_clusters)
scRNASeJeAP@meta.data$cellType2<- paste0(scRNASeJeAP@meta.data$cellType,scRNASeJeAP@meta.data$seurat_clusters)
library(dplyr)
 
# 将df1转换为tibble，以便使用left_join，并设置行名为一个临时列
scRNASeJeAP_tibble <- tibble::rownames_to_column(as_tibble(scRNASeJeAP@meta.data), var = "rowname")
scRNASeJeAP_tibble$rowname=rownames(scRNASeJeAP@meta.data)
SeuratObject_tibble <- tibble::rownames_to_column(as_tibble(SeuratObject@meta.data), var = "rowname")
SeuratObject_tibble$rowname=rownames(SeuratObject@meta.data)
# 使用left_join合并数据框，然后根据结果替换df2中的celltype2列
merged_df <- left_join(SeuratObject_tibble, scRNASeJeAP_tibble, by = "rowname")
 
# 检查是否有NA值，如果有，则用原始的celltype2值填充
merged_df$cellType2 <- as.character(merged_df$cellType2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
# 使用ifelse进行替换
merged_df$celltype2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)
 
# 如果需要，将celltype2转换回因子（并更新水平集）
# 注意：这一步只在celltype2原本就是因子且你需要保持它为因子类型时才需要
# levels_to_use <- unique(c(as.character(merged_df$celltype2), levels(as.factor(scRNASeJeAP@meta.data$cellType2))))
merged_df$celltype2 <- factor(merged_df$celltype2, levels = c('AP0','AP1','AP2','AP3','AP4','AP5','AP6','AP7','myoFB'))
SeuratObject@meta.data$celltype2=merged_df$celltype2

SeuratObject2=SeuratObject
DefaultAssay(SeuratObject2) <- "RNA"

SeuratObject2 <- NormalizeData(SeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

SeuratObject2 <- FindVariableFeatures(SeuratObject2,nfeatures = 3000)
SeuratObject2 <- ScaleData(SeuratObject2,vars.to.regress = 'percent.mt')
SeuratObject2 <- RunPCA(SeuratObject2)
SeuratObject2 <- FindNeighbors(SeuratObject2, reduction = "pca", dims = 1:30)
SeuratObject2 <- FindClusters(SeuratObject2, resolution = 0.3)
SeuratObject2 <- RunUMAP(SeuratObject2, reduction = "pca", dims = 1:30) 
library(harmony)
SeuratObject2 <- RunHarmony(SeuratObject2, group.by.vars = "orig.ident")
SeuratObject2  <- FindNeighbors(SeuratObject2 , dims = 1:30 , reduction = "harmony")
SeuratObject2  <- FindClusters(SeuratObject2, save.snn=T , resolution = 0.3)
SeuratObject2  <- RunUMAP(SeuratObject2 , dims=1:30,reduction='harmony')

p1=DimPlot(SeuratObject2,group.by = "orig.ident")
# 查看聚类信息
p2=DimPlot(SeuratObject2,group.by = "seurat_clusters",label=T)
p3=DimPlot(SeuratObject2,group.by = "celltype2",label=T)
combine<-CombinePlots(list(p1,p2))
ggsave(filename = "sejep_DimPlot_sample_seurat_clusters.pdf", plot = combine, height = 6, width = 10)

library(monocle3)
SeuratObject=SeuratObject2
data <- GetAssayData(SeuratObject, assay = 'RNA', slot = 'counts')
cell_metadata <- SeuratObject@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 19)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP",color_cells_by="celltype2") + ggtitle('cds.umap')
########可视化部分基因
# ciliated_genes <- c("MITF",
#                     "SOX10",
#                     "AXL",
#                     "NGFR")

# plot_cells(cds,
#            genes=ciliated_genes,
#            label_cell_groups=FALSE,
#            show_trajectory_graph=FALSE)
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(SeuratObject, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype2") + ggtitle('int.umap')
# p1+p2
combine<-CombinePlots(list(p1,p2))
ggsave(filename = "sejecds_int.umap.pdf", plot = combine, height = 6, width = 8)

expr_matrix <- as.data.frame(SeuratObject@assays$RNA@counts)
# 细胞注释矩阵（列为细胞名）
sample_sheet <- SeuratObject@meta.data
# 基因注释矩阵（行为基因名）
gene_annotation <- data.frame(
  gene_short_name=rownames(SeuratObject@assays$RNA),
  row.names = rownames(SeuratObject@assays$RNA)
)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

#创建monocle对象
cd <- newCellDataSet(as(as.matrix(expr_matrix),"sparseMatrix"), 
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

###计算数据的size factors和dispersions
cd <- estimateSizeFactors(cd)
cd <- estimateDispersions(cd)
# 过滤低于1%细胞中检出的基因，最低表达阈值为0.5
cd <- detectGenes(cd, min_expr = 0.5)
expressed_genes <- row.names(subset(fData(cd), num_cells_expressed > nrow(sample_sheet) * 0.01))
# 查看筛选后的基因个数（用于后面基因的筛选）
length(expressed_genes)

###构建轨迹
disp_table <- dispersionTable(cd)
# 高度离散基因的筛选标准，可根据数据情况设置mean_expression的值
ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)
# 根据数据的某种分类进行差异分析
diff_test_res <- differentialGeneTest(cd[expressed_genes,],fullModelFormulaStr = "~celltype2")
# 筛选差异基因（q < 1e-5并且属于之前计算的expressed_genes列表中）
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-5))
ordering_genes <- intersect(ordering_genes, expressed_genes)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)
# 选取的基因数目为每个细胞的维度，基于默认的'DDRTree'方法进行数据降维
cd <- reduceDimension(cd, max_components = 2, method = 'DDRTree')
# 对细胞进行排序，由于排序无法区分起点和终点，若分析所得时序与实际相反，根据“reverse”参数进行调整，默认reverse=F
cd <- orderCells(cd, reverse = F) 
# 排序好的细胞可以进行可视化，可标注细胞的各注释信息)
# 查看细胞注释信息
head(cd@phenoData@data)
saveRDS(cd,file='SeJem2_cd_apmyofb.rds')
cd <- readRDS('SeJem2_cd_apmyofb.rds')
# 设置颜色
color1 <- c(brewer.pal(8, "Set1"))
getPalette <- colorRampPalette(brewer.pal(6, "Set1"))
# # 不拆分
# # 按照表型进行映射
# pdf(file="cd_cellType2.trajectory.pdf",width=6.5,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "celltype2") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"celltype2"]))))
# dev.off()
# # 按照表型进行映射
# pdf(file="cd_cellType.trajectory.pdf",width=6.5,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"celltype2"]))))
# dev.off()

# 拆分
pData(cd)$sample=SeuratObject$orig.ident
nrows_state=1
pdf(file="SeJecd_cellType.trajectory.pdf",width=6,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"cellType"]))))
plot_cell_trajectory(cd, color_by = "cellType", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~sample, nrow = nrows_state, scales = "free")+ scale_color_manual(values=c('red','blue'))+ggtitle("cellType") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()

0,2,6,8,1,3,4,5,7,myoFB
pData(cd)$celltype2 <- factor(pData(cd)$celltype2, levels = c('AP0','AP2','AP6','AP8','AP1','AP3','AP4','AP5','AP7','myoFB'))
pdf(file="SeJecd_sample.trajectory.pdf",width=6,height=12)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"cellType"]))))
plot_cell_trajectory(cd, color_by = "sample", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~celltype2, nrow = 2, scales = "free")+ggtitle("sample") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()


scRNASeJeAP=readRDS('../cellchat/scRNASeJeAPharmony.rda')
Idents(scRNASeJeAP) <- "seurat_clusters"
markers <- FindMarkers(scRNASeJeAP, ident.1 = 6, only.pos = TRUE)
markers %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 8) %>%
    ungroup() -> top8
write.csv(file='AP6.markers.csv',markers)
DefaultAssay(scRNASeJeAP) <- "RNA"
pdf(file="sejeap6FeaturePlottop10gene.pdf",,height = 6, width = 6)
FeaturePlot(scRNASeJeAP, features = c(rownames(top8),'Nucb2'))
dev.off()

pdf(file="sejeap6dotplotaddNucb2.pdf",,height = 6, width = 6)
DotPlot(scRNASeJeAP, features = c(rownames(top8),'Nucb2'))+ coord_flip()
dev.off()


scRNA=readRDS('../scRNA_seurat_cellType.RDS')
scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c("Je","Ji","JT","Se","Si","ST"))
scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))
scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
scRNASTJT=subset(x = scRNA, subset = (orig.ident == "ST"|orig.ident == "JT"))
samples=c('SeJe','SiJi','STJT')
features<-c('Nucb2')
DefaultAssay(scRNASeJe) <- "RNA"
DefaultAssay(scRNASiJi) <- "RNA"
DefaultAssay(scRNASTJT) <- "RNA"
i=1
for (obj in c(scRNASeJe,scRNASiJi,scRNASTJT)){
  sample=samples[i]
  pdf(file=sprintf("Nucb2_%s.pdf",sample),height = 6, width = 9)
  par(mfrow = c(2,1), xpd=TRUE)
  # FeaturePlot(object = obj, features = features,split.by = 'orig.ident')
  # print(FeaturePlot(object = obj, features = features,split.by = 'orig.ident'))
  DotPlot(obj,group.by = 'cellType', features = unique(features),split.by = 'orig.ident') + coord_flip()+RotatedAxis()+ theme(text = element_text(size = 30))+theme(axis.text.x = element_text(face="bold", size=9),
           axis.text.y = element_text(face="bold",  size=20))
  print(DotPlot(obj,group.by = 'cellType', features = unique(features),split.by = 'orig.ident') + coord_flip()+RotatedAxis()+ theme(text = element_text(size = 30))+theme(axis.text.x = element_text(face="bold", size=9),
           axis.text.y = element_text(face="bold",  size=20)))
  dev.off()
  i=i+1
}


scRNASeJeAP=readRDS('../cellchat/scRNASeJeAPharmony.rda')

scRNASeJeAP@meta.data$cellType2<- 'APNo6'
Idents(scRNASeJeAP)='seurat_clusters'
filtered_cellsAP6 <- WhichCells(scRNASeJeAP,idents='6')
Idents(scRNASeJeAP)='cellType2'
Idents(scRNASeJeAP, cells = filtered_cellsAP6) <- 'AP6'
table(Idents(scRNASeJeAP))
scRNASeJeAP@meta.data$cellType2=Idents(scRNASeJeAP)
library(limma)
library(edgeR)
data=scRNASeJeAP@assays$RNA@data
dd_CM<- as.data.frame(as.matrix(data))
# colnames(dd_CM)

group2<-  factor(scRNASeJeAP@meta.data$cellType2,levels = c("APNo6", "AP6"))#顺序靠后的Se vs 顺序靠前的Je;#对照组在前，处理组在后
table(group2)
dgelist <- DGEList(counts =dd_CM, group = group2)
#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
#（3）标准化，以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#差异表达基因分析
design <- model.matrix(~group2)
#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
sample1='AP6'
sample2='APNo6'
write.table(lrt, sprintf('%svs%s.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)
library(ggplot2)

d1<- read.table(sprintf('%svs%s.txt',sample1,sample2))
colnames(d1)
d1$gene<- rownames(d1)

# 设置p_value和logFC的阈值
cut_off_FDR = 0.05  #统计显著性
cut_off_logFC = 0.5           #差异倍数值
d1$change<- 'Stable'
# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
d1$change = ifelse(d1$FDR< cut_off_FDR & abs(d1$logFC) > cut_off_logFC, 
                        ifelse(d1$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(d1)
table(d1$change)

p <- ggplot(
  # 数据、映射、颜色
  d1, aes(x = logFC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_FDR),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  scale_x_continuous(limits = c(-5, 5))+
  scale_y_continuous(limits = c(0,30))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('AP6 v.s. APNo6')
# p
ggsave(filename = "diffgene2nolabel.pdf", height = 6, width = 8, plot = p)
# 将需要标记的基因放置在label列(logFC >= 5)
library(ggrepel)
d1$label <- ifelse(d1$PValue < cut_off_FDR & abs(d1$logFC) >= 3,
                        as.character(d1$gene), "")
d1$label[d1$gene=='Nucb2']='Nucb2'

p=p + geom_label_repel(data = d1, aes(x = d1$logFC, 
                                         y = -log10(d1$PValue), 
                                         label = label),
                     size = 3, box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"), 
                     segment.color = "black", 
                     show.legend = FALSE,max.overlaps = Inf)
ggsave(filename = "diffgene2.pdf", height = 6, width = 6, plot = p)


# 
scRNASiJiAP=readRDS('../APSIJI/scRNASiJiAPharmony.rda')

scRNASiJiAP@meta.data$cellType2<- 'APNo1'
Idents(scRNASiJiAP)='seurat_clusters'
filtered_cellsAP6 <- WhichCells(scRNASiJiAP,idents='1')
Idents(scRNASiJiAP)='cellType2'
Idents(scRNASiJiAP, cells = filtered_cellsAP6) <- 'AP1'
table(Idents(scRNASiJiAP))
scRNASiJiAP@meta.data$cellType2=Idents(scRNASiJiAP)
library(limma)
library(edgeR)
data=scRNASiJiAP@assays$RNA@data
dd_CM<- as.data.frame(as.matrix(data))
# colnames(dd_CM)

group2<-  factor(scRNASiJiAP@meta.data$cellType2,levels = c("APNo1", "AP1"))#顺序靠后的Se vs 顺序靠前的Je;#对照组在前，处理组在后
table(group2)
dgelist <- DGEList(counts =dd_CM, group = group2)
#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
#（3）标准化，以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#差异表达基因分析
design <- model.matrix(~group2)
#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
sample1='AP1'
sample2='APNo1'
write.table(lrt, sprintf('%svs%s.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)
library(ggplot2)

d1<- read.table(sprintf('%svs%s.txt',sample1,sample2))
colnames(d1)
d1$gene<- rownames(d1)

# 设置p_value和logFC的阈值
cut_off_FDR = 0.05  #统计显著性
cut_off_logFC = 0.5           #差异倍数值
d1$change<- 'Stable'
# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
d1$change = ifelse(d1$FDR< cut_off_FDR & abs(d1$logFC) > cut_off_logFC, 
                        ifelse(d1$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(d1)
table(d1$change)
write.csv(file='sijiAP1diff.csv',d1)
p <- ggplot(
  # 数据、映射、颜色
  d1, aes(x = logFC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_FDR),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  scale_x_continuous(limits = c(-2.5, 2.5))+
  scale_y_continuous(limits = c(0,10))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('AP1 v.s. APNo1')
# p
ggsave(filename = "sijidiffgene2nolabelAP1.pdf", height = 6, width = 6, plot = p)
# 将需要标记的基因放置在label列(logFC >= 5)
library(ggrepel)
d1$label <- ifelse(d1$PValue < cut_off_FDR & abs(d1$logFC) >= 3,
                        as.character(d1$gene), "")
# d1$label[d1$gene=='Nucb2']='Nucb2'

p=p + geom_label_repel(data = d1, aes(x = d1$logFC, 
                                         y = -log10(d1$PValue), 
                                         label = label),
                     size = 3, box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"), 
                     segment.color = "black", 
                     show.legend = FALSE,max.overlaps = Inf)
ggsave(filename = "sijidiffgene2AP1.pdf", height = 6, width = 6, plot = p)

genes=rownames(d1)[d1$change=='Down']

scRNA <- readRDS("../scRNA_seurat_cellType.RDS")
scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
scRNASiJiAPmyoFBEC=subset(x = scRNASiJi, subset = (cellType == "AP"|cellType == "myoFB"|cellType == "EC"))
SeuratObject=scRNASiJiAPmyoFBEC

scRNASiJiAP=readRDS('../APSIJI/scRNASiJiAPharmony.rda')
scRNASiJiAP@meta.data$cellType2<- 'APNo1'
Idents(scRNASiJiAP)='seurat_clusters'
filtered_cellsAP6 <- WhichCells(scRNASiJiAP,idents='1')
Idents(scRNASiJiAP)='cellType2'
Idents(scRNASiJiAP, cells = filtered_cellsAP6) <- 'AP1'
table(Idents(scRNASiJiAP))
scRNASiJiAP@meta.data$cellType2=Idents(scRNASiJiAP)

library(dplyr)
 
# 将df1转换为tibble，以便使用left_join，并设置行名为一个临时列
scRNASiJiAP_tibble <- tibble::rownames_to_column(as_tibble(scRNASiJiAP@meta.data), var = "rowname")
scRNASiJiAP_tibble$rowname=rownames(scRNASiJiAP@meta.data)
SeuratObject_tibble <- tibble::rownames_to_column(as_tibble(SeuratObject@meta.data), var = "rowname")
SeuratObject_tibble$rowname=rownames(SeuratObject@meta.data)
# 使用left_join合并数据框，然后根据结果替换df2中的celltype2列
merged_df <- left_join(SeuratObject_tibble, scRNASiJiAP_tibble, by = "rowname")
 
# 检查是否有NA值，如果有，则用原始的celltype2值填充
merged_df$cellType2 <- as.character(merged_df$cellType2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
# 使用ifelse进行替换
merged_df$celltype2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)
 
# 如果需要，将celltype2转换回因子（并更新水平集）
# 注意：这一步只在celltype2原本就是因子且你需要保持它为因子类型时才需要
# levels_to_use <- unique(c(as.character(merged_df$celltype2), levels(as.factor(scRNASeJeAP@meta.data$cellType2))))
merged_df$celltype2 <- factor(merged_df$celltype2, levels = c('AP1','APNo1','myoFB','EC'))
SeuratObject@meta.data$celltype2=merged_df$celltype2

dot=DotPlot(SeuratObject, features = c(genes),group.by='celltype2')
dotdata=dot$data
write.csv(file='sijidotdata.csv',dotdata)
write.csv(file='sijidotdataAPNo1.csv',dotdata[dotdata$id=='APNo1',])
write.csv(file='sijidotdataAP1.csv',dotdata[dotdata$id=='AP1',])
write.csv(file='sijidotdatamyoFB.csv',dotdata[dotdata$id=='myoFB',])
write.csv(file='sijidotdataEC.csv',dotdata[dotdata$id=='EC',])



scRNA <- readRDS("../scRNA_seurat_cellType.RDS")
scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
scRNASiJiAP=subset(x = scRNASiJi, subset = (cellType == "AP"))
SeuratObject=scRNASiJiAP
DefaultAssay(SeuratObject) <- "RNA"

SeuratObject2 <- NormalizeData(SeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

SeuratObject2 <- FindVariableFeatures(SeuratObject2,nfeatures = 3000)
SeuratObject2 <- ScaleData(SeuratObject2,vars.to.regress = 'percent.mt')
SeuratObject2 <- RunPCA(SeuratObject2)
SeuratObject2 <- FindNeighbors(SeuratObject2, reduction = "pca", dims = 1:30)
SeuratObject2 <- FindClusters(SeuratObject2, resolution = 0.3)
SeuratObject2 <- RunUMAP(SeuratObject2, reduction = "pca", dims = 1:30) 
library(harmony)
SeuratObject2 <- RunHarmony(SeuratObject2, group.by.vars = "orig.ident")
SeuratObject2  <- FindNeighbors(SeuratObject2 , dims = 1:30 , reduction = "harmony")
SeuratObject2  <- FindClusters(SeuratObject2, save.snn=T , resolution = 0.3)
SeuratObject2  <- RunUMAP(SeuratObject2 , dims=1:30,reduction='harmony')

p1=DimPlot(SeuratObject2,group.by = "orig.ident")
# 查看聚类信息
p2=DimPlot(SeuratObject2,group.by = "seurat_clusters",label=T)
# p3=DimPlot(SeuratObject2,group.by = "celltype2",label=T)
combine<-CombinePlots(list(p1,p2))
ggsave(filename = "sijip_DimPlot_sample_seurat_clusters.pdf", plot = combine, height = 6, width = 10)
saveRDS(SeuratObject2,file= "sijiap0.3.rds")


scRNA=readRDS('../scRNA_seurat_cellType.RDS')
scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
scRNASiJiAPmyoFB=subset(x=scRNASiJi,subset=(cellType == "AP"|cellType == "myoFB"))

SeuratObject=scRNASiJiAPmyoFB
scRNASiJiAP=readRDS('sijiap0.2.rds')
scRNASiJiAP$seurat_clusters=Idents(scRNASiJiAP)
# object=scRNASiJiAP
# object$seurat_clusters <- as.character(object$seurat_clusters)
# object$seurat_clusters[object$seurat_clusters == "5"] <- "1"
# object$seurat_clusters[object$seurat_clusters == "8"] <- "5"
# object$seurat_clusters <- factor(object$seurat_clusters)

# SeuratObject$Cluster<- paste0('cluster',SeuratObject$seurat_clusters)
scRNASiJiAP@meta.data$cellType2<- paste0(scRNASiJiAP@meta.data$cellType,scRNASiJiAP@meta.data$seurat_clusters)
library(dplyr)
 
# 将df1转换为tibble，以便使用left_join，并设置行名为一个临时列
scRNASiJiAP_tibble <- tibble::rownames_to_column(as_tibble(scRNASiJiAP@meta.data), var = "rowname")
scRNASiJiAP_tibble$rowname=rownames(scRNASiJiAP@meta.data)
SeuratObject_tibble <- tibble::rownames_to_column(as_tibble(SeuratObject@meta.data), var = "rowname")
SeuratObject_tibble$rowname=rownames(SeuratObject@meta.data)
# 使用left_join合并数据框，然后根据结果替换df2中的celltype2列
merged_df <- left_join(SeuratObject_tibble, scRNASiJiAP_tibble, by = "rowname")
 
# 检查是否有NA值，如果有，则用原始的celltype2值填充
merged_df$cellType2 <- as.character(merged_df$cellType2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
# 使用ifelse进行替换
merged_df$celltype2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)
 
# 如果需要，将celltype2转换回因子（并更新水平集）
# 注意：这一步只在celltype2原本就是因子且你需要保持它为因子类型时才需要
levels_to_use <- unique(c(as.character(merged_df$celltype2), levels(as.factor(scRNASiJiAP@meta.data$cellType))))
merged_df$celltype2 <- factor(merged_df$celltype2, levels = levels_to_use)
SeuratObject@meta.data$celltype2=merged_df$celltype2

SeuratObject2=SeuratObject
DefaultAssay(SeuratObject2) <- "RNA"

SeuratObject2 <- NormalizeData(SeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

SeuratObject2 <- FindVariableFeatures(SeuratObject2,nfeatures = 3000)
SeuratObject2 <- ScaleData(SeuratObject2,vars.to.regress = 'percent.mt')
SeuratObject2 <- RunPCA(SeuratObject2)
SeuratObject2 <- FindNeighbors(SeuratObject2, reduction = "pca", dims = 1:30)
SeuratObject2 <- FindClusters(SeuratObject2, resolution = 0.2)
SeuratObject2 <- RunUMAP(SeuratObject2, reduction = "pca", dims = 1:30) 
library(harmony)
SeuratObject2 <- RunHarmony(SeuratObject2, group.by.vars = "orig.ident")
SeuratObject2  <- FindNeighbors(SeuratObject2 , dims = 1:30 , reduction = "harmony")
SeuratObject2  <- FindClusters(SeuratObject2, save.snn=T , resolution = 0.2)
SeuratObject2  <- RunUMAP(SeuratObject2 , dims=1:30,reduction='harmony')

p1=DimPlot(SeuratObject2,group.by = "orig.ident")
# 查看聚类信息
p2=DimPlot(SeuratObject2,group.by = "seurat_clusters",label=T)
p3=DimPlot(SeuratObject2,group.by = "celltype2",label=T)
combine<-CombinePlots(list(p1,p2,p3))
ggsave(filename = "sijip_DimPlot_sample_seurat_clusters3.pdf", plot = combine, height = 6, width = 15)


library(monocle)
detach("package:monocle3", unload = TRUE)

scRNA <- readRDS("../scRNA_seurat_cellType.RDS")
scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
scRNASiJiAPmyoFB=subset(x = scRNASiJi, subset = (cellType == "AP"|cellType == "myoFB"))
SeuratObject=scRNASiJiAPmyoFB
scRNASiJiAP=readRDS('sijiap0.2.rds')
scRNASiJiAP$seurat_clusters=Idents(scRNASiJiAP)
# SeuratObject$Cluster<- paste0('cluster',SeuratObject$seurat_clusters)
scRNASiJiAP@meta.data$cellType2<- paste0(scRNASiJiAP@meta.data$cellType,scRNASiJiAP@meta.data$seurat_clusters)
library(dplyr)
 
# 将df1转换为tibble，以便使用left_join，并设置行名为一个临时列
scRNASiJiAP_tibble <- tibble::rownames_to_column(as_tibble(scRNASiJiAP@meta.data), var = "rowname")
scRNASiJiAP_tibble$rowname=rownames(scRNASiJiAP@meta.data)
SeuratObject_tibble <- tibble::rownames_to_column(as_tibble(SeuratObject@meta.data), var = "rowname")
SeuratObject_tibble$rowname=rownames(SeuratObject@meta.data)
# 使用left_join合并数据框，然后根据结果替换df2中的celltype2列
merged_df <- left_join(SeuratObject_tibble, scRNASiJiAP_tibble, by = "rowname")
 
# 检查是否有NA值，如果有，则用原始的celltype2值填充
merged_df$cellType2 <- as.character(merged_df$cellType2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
# 使用ifelse进行替换
merged_df$celltype2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)
 
# 如果需要，将celltype2转换回因子（并更新水平集）
# 注意：这一步只在celltype2原本就是因子且你需要保持它为因子类型时才需要
merged_df$celltype2 <- factor(merged_df$celltype2, levels = c('AP0','AP1','AP2','AP3','AP4','AP5','AP6','myoFB'))
SeuratObject@meta.data$celltype2=merged_df$celltype2
expr_matrix <- as.data.frame(SeuratObject@assays$RNA@counts)
# 细胞注释矩阵（列为细胞名）
sample_sheet <- SeuratObject@meta.data
# 基因注释矩阵（行为基因名）
gene_annotation <- data.frame(
  gene_short_name=rownames(SeuratObject@assays$RNA),
  row.names = rownames(SeuratObject@assays$RNA)
)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

#创建monocle对象
cd <- newCellDataSet(as(as.matrix(expr_matrix),"sparseMatrix"), 
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

###计算数据的size factors和dispersions
cd <- estimateSizeFactors(cd)
cd <- estimateDispersions(cd)
# 过滤低于1%细胞中检出的基因，最低表达阈值为0.5
cd <- detectGenes(cd, min_expr = 0.5)
expressed_genes <- row.names(subset(fData(cd), num_cells_expressed > nrow(sample_sheet) * 0.01))
# 查看筛选后的基因个数（用于后面基因的筛选）
length(expressed_genes)

###构建轨迹
disp_table <- dispersionTable(cd)
# 高度离散基因的筛选标准，可根据数据情况设置mean_expression的值
ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)
# 根据数据的某种分类进行差异分析
diff_test_res <- differentialGeneTest(cd[expressed_genes,],fullModelFormulaStr = "~celltype2")
# 筛选差异基因（q < 1e-5并且属于之前计算的expressed_genes列表中）
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-5))
ordering_genes <- intersect(ordering_genes, expressed_genes)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)
# 选取的基因数目为每个细胞的维度，基于默认的'DDRTree'方法进行数据降维
cd <- reduceDimension(cd, max_components = 2, method = 'DDRTree')
# 对细胞进行排序，由于排序无法区分起点和终点，若分析所得时序与实际相反，根据“reverse”参数进行调整，默认reverse=F
cd <- orderCells(cd, reverse = F) 
# 排序好的细胞可以进行可视化，可标注细胞的各注释信息)
# 查看细胞注释信息
head(cd@phenoData@data)
saveRDS(cd,file='SiJim2_cd_apmyofb.rds')
cd <- readRDS('SiJim2_cd_apmyofb.rds')
# 设置颜色
color1 <- c(brewer.pal(8, "Set1"))
getPalette <- colorRampPalette(brewer.pal(6, "Set1"))
# # 不拆分
# # 按照表型进行映射
# pdf(file="cd_cellType2.trajectory.pdf",width=6.5,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "celltype2") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"celltype2"]))))
# dev.off()
# # 按照表型进行映射
# pdf(file="cd_cellType.trajectory.pdf",width=6.5,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"celltype2"]))))
# dev.off()

cd=readRDS('SiJim2_cd_apmyofb.rds')
# 拆分
pData(cd)$sample=SeuratObject$orig.ident
pData(cd)$sample=factor(pData(cd)$sample, levels = c('Ji','Si'))
nrows_state=1
pdf(file="SiJicd_cellType.trajectory.pdf",width=12,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"cellType"]))))
plot_cell_trajectory(cd, color_by = "cellType", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~sample, nrow = nrows_state, scales = "free")+ scale_color_manual(values=c('red','blue'))+ggtitle("cellType") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()

0,2,6,8,1,3,4,5,7,myoFB
pData(cd)$celltype2 <- factor(pData(cd)$celltype2, levels = c('AP0','AP1','AP6','AP5','AP3','AP4','AP2','myoFB'))
pdf(file="SiJicd_sample.trajectory.pdf",width=12,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"cellType"]))))
plot_cell_trajectory(cd, color_by = "sample", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~celltype2, nrow = 3, scales = "free")+ggtitle("sample") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()


# 
scRNASiJiAP=readRDS('sijiap0.2.rds')

scRNASiJiAP@meta.data$cellType2<- 'APNo2'
Idents(scRNASiJiAP)='seurat_clusters'
filtered_cellsAP6 <- WhichCells(scRNASiJiAP,idents='2')
Idents(scRNASiJiAP)='cellType2'
Idents(scRNASiJiAP, cells = filtered_cellsAP6) <- 'AP2'
table(Idents(scRNASiJiAP))
scRNASiJiAP@meta.data$cellType2=Idents(scRNASiJiAP)
library(limma)
library(edgeR)
data=scRNASiJiAP@assays$RNA@data
dd_CM<- as.data.frame(as.matrix(data))
# colnames(dd_CM)

group2<-  factor(scRNASiJiAP@meta.data$cellType2,levels = c("APNo2", "AP2"))#顺序靠后的Se vs 顺序靠前的Je;#对照组在前，处理组在后
table(group2)
dgelist <- DGEList(counts =dd_CM, group = group2)
#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
#（3）标准化，以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#差异表达基因分析
design <- model.matrix(~group2)
#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
sample1='AP2'
sample2='APNo2'
write.table(lrt, sprintf('%svs%s.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)
library(ggplot2)

d1<- read.table(sprintf('%svs%s.txt',sample1,sample2))
colnames(d1)
d1$gene<- rownames(d1)

# 设置p_value和logFC的阈值
cut_off_FDR = 0.05  #统计显著性
cut_off_logFC = 0.5           #差异倍数值
d1$change<- 'Stable'
# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
d1$change = ifelse(d1$FDR< cut_off_FDR & abs(d1$logFC) > cut_off_logFC, 
                        ifelse(d1$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(d1)
table(d1$change)
write.csv(file='sijiAP2diff.csv',d1)
p <- ggplot(
  # 数据、映射、颜色
  d1, aes(x = logFC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_FDR),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  scale_x_continuous(limits = c(-2.5, 2.5))+
  scale_y_continuous(limits = c(0,10))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('AP2 v.s. APNo2')
# p
ggsave(filename = "sijidiffgene2nolabelAP2.pdf", height = 6, width = 6, plot = p)
# 将需要标记的基因放置在label列(logFC >= 5)
library(ggrepel)
d1$label <- ifelse(d1$PValue < cut_off_FDR & abs(d1$logFC) >= 3,
                        as.character(d1$gene), "")
# d1$label[d1$gene=='Nucb2']='Nucb2'

p=p + geom_label_repel(data = d1, aes(x = d1$logFC, 
                                         y = -log10(d1$PValue), 
                                         label = label),
                     size = 3, box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"), 
                     segment.color = "black", 
                     show.legend = FALSE,max.overlaps = Inf)
ggsave(filename = "sijidiffgene2AP2.pdf", height = 6, width = 6, plot = p)

genes=rownames(d1)[d1$change=='Up']
scRNA <- readRDS("../scRNA_seurat_cellType.RDS")
scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
scRNASiJiAPmyoFBEC=subset(x = scRNASiJi, subset = (cellType == "AP"|cellType == "myoFB"|cellType == "EC"))
SeuratObject=scRNASiJiAPmyoFBEC

scRNASiJiAP=readRDS('')
scRNASiJiAP@meta.data$cellType2<- 'APNo2'
Idents(scRNASiJiAP)='seurat_clusters'
filtered_cellsAP6 <- WhichCells(scRNASiJiAP,idents='2')
Idents(scRNASiJiAP)='cellType2'
Idents(scRNASiJiAP, cells = filtered_cellsAP6) <- 'AP2'
table(Idents(scRNASiJiAP))
scRNASiJiAP@meta.data$cellType2=Idents(scRNASiJiAP)

library(dplyr)
 
# 将df1转换为tibble，以便使用left_join，并设置行名为一个临时列
scRNASiJiAP_tibble <- tibble::rownames_to_column(as_tibble(scRNASiJiAP@meta.data), var = "rowname")
scRNASiJiAP_tibble$rowname=rownames(scRNASiJiAP@meta.data)
SeuratObject_tibble <- tibble::rownames_to_column(as_tibble(SeuratObject@meta.data), var = "rowname")
SeuratObject_tibble$rowname=rownames(SeuratObject@meta.data)
# 使用left_join合并数据框，然后根据结果替换df2中的celltype2列
merged_df <- left_join(SeuratObject_tibble, scRNASiJiAP_tibble, by = "rowname")
 
# 检查是否有NA值，如果有，则用原始的celltype2值填充
merged_df$cellType2 <- as.character(merged_df$cellType2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
# 使用ifelse进行替换
merged_df$celltype2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)
 
# 如果需要，将celltype2转换回因子（并更新水平集）
# 注意：这一步只在celltype2原本就是因子且你需要保持它为因子类型时才需要
# levels_to_use <- unique(c(as.character(merged_df$celltype2), levels(as.factor(scRNASeJeAP@meta.data$cellType2))))
merged_df$celltype2 <- factor(merged_df$celltype2, levels = c('AP2','APNo2','myoFB','EC'))
SeuratObject@meta.data$celltype2=merged_df$celltype2

dot=DotPlot(SeuratObject, features = c(genes),group.by='celltype2')
dotdata=dot$data
write.csv(file='sijidotdata.csv',dotdata)
write.csv(file='sijidotdataAPNo5.csv',dotdata[dotdata$id=='APNo2',])
write.csv(file='sijidotdataAP5.csv',dotdata[dotdata$id=='AP2',])
write.csv(file='sijidotdatamyoFB.csv',dotdata[dotdata$id=='myoFB',])
write.csv(file='sijidotdataEC.csv',dotdata[dotdata$id=='EC',])


degs=markers
degs_down <- degs[degs$avg_log2FC < -0.25 & degs$p_val_adj < 0.05, ]
degs_down$gene=rownames(degs_down)
library(org.Mm.eg.db)
enrich_go <- enrichGO(gene = degs_down$gene, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")
p=dotplot(enrich_go)
ggsave(filename = "sijienrich.pdf", height = 6, width = 6, plot = p)
p=FeaturePlot(scRNASiJiAP, features = c("HSPA6", "DDIT4","BBC3"), reduction = "umap")
ggsave(filename = "sijienrichfeature.pdf", height = 6, width = 6, plot = p)

scRNASiJiAP=readRDS('../APSIJI/scRNASiJiAPharmony.rda')
markers <- FindMarkers(scRNASiJiAP, ident.1 = 1,logfc.threshold=0)
d1=markers


d1$gene<- rownames(d1)
cut_off_FDR = 0.05  #统计显著性
cut_off_logFC = 0.5           #差异倍数值
d1$change<- 'Stable'
d1$change = ifelse(d1$p_val_adj< cut_off_FDR & abs(d1$avg_log2FC) > cut_off_logFC, 
                        ifelse(d1$avg_log2FC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(d1)
table(d1$change)

p <- ggplot(
  # 数据、映射、颜色
  d1, aes(x = avg_log2FC, y = -log10(p_val_adj), colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_FDR),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  # scale_x_continuous(limits = c(-2.5, 2.5))+
  # scale_y_continuous(limits = c(0,10))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('AP1 v.s. APNo1')
# p
ggsave(filename = "sijidiffgene2nolabelAP1markers.pdf", height = 6, width = 6, plot = p)



scRNASeJeAP=readRDS('../cellchat/scRNASeJeAPharmony.rda')
markers2 <- FindMarkers(scRNASeJeAP, ident.1 = 6,logfc.threshold=0)
d1=markers2
d1$gene<- rownames(d1)
cut_off_FDR = 0.05  #统计显著性
cut_off_logFC = 0.5           #差异倍数值
d1$change<- 'Stable'
d1$change = ifelse(d1$p_val< cut_off_FDR & abs(d1$avg_log2FC) > cut_off_logFC, 
                        ifelse(d1$avg_log2FC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(d1)
table(d1$change)

p <- ggplot(
  # 数据、映射、颜色
  d1, aes(x = avg_log2FC, y = -log10(p_val), colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_FDR),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  scale_x_continuous(limits = c(-5, 5))+
  scale_y_continuous(limits = c(0,30))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('AP6 v.s. APNo6')
# p
ggsave(filename = "sejediffgene2nolabelAP6markersrange.pdf", height = 6, width = 6, plot = p)