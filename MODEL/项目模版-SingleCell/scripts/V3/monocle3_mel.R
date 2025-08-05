############创建seurat对象
library(Seurat)
library(monocle3)
# devtools::install_github('cole-trapnell-lab/monocle3')
# devtools::install_github('NameAIAK/Cellchat-works')
library(tidyverse)
library(patchwork)
# library(harmony)
# BiocManager::install("harmony")
##运行mel.R

##########################monocle3轨迹分析
##创建CDS对象并预处理数据
scRNA <- readRDS("../scRNA_seurat_cellType.RDS")
scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))
# scRNASeJeAPmyoFBECFB=subset(x = scRNASeJe, subset = (cellType == "AP"|cellType == "myoFB"|cellType == "EC"|cellType == "FB"))
# SeuratObject=scRNASeJeAPmyoFBECFB
scRNASeJeAPmyoFB=subset(x = scRNASeJe, subset = (cellType == "AP"))
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
ggsave(filename = "p_DimPlot_sample_seurat_clusters.pdf", plot = combine, height = 6, width = 7)


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
ggsave(filename = "cds_int.umap.pdf", plot = combine, height = 6, width = 7)
## Monocle3聚类分区
cds <- cluster_cells(cds,resolution = 0.1)
p1 <- plot_cells(cds, color_cells_by='celltype2',show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
# p = wrap_plots(p1, p2)

## 识别轨迹
cds <- learn_graph(cds)
p3 = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = T, color_cells_by='celltype2', 
               label_branch_points = T)


##细胞按拟时排序
##############################################################
# cds <- order_cells(cds) 存在bug，使用辅助线选择root细胞
# options(browser="firefox")
# cds <- order_cells(cds)
# p4=plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
#            label_leaves = FALSE,  label_branch_points = FALSE)
combine<-CombinePlots(list(p1,p2,p3))
ggsave(filename = "cds_int.umap2.pdf", plot = combine, height = 6, width = 7)
#############################################################
##############monocle3差异分析
##寻找拟时轨迹差异基因
#graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
#空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=3)
#挑选top10画图展示
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#基因表达趋势图
plot_genes_in_pseudotime(cds[Track_genes_sig,], 
                         min_expr=0.5, ncol = 2)
#FeaturePlot图
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
##寻找共表达模块
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)$predicted.id)
cell_group[,2]=gatk[,1]
agg_mat <- aggregate_gene_expression(cds,gene_group_df=gene_module, cell_group_df=cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")



############创建seurat对象
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(harmony)

##运行mel.R

##########################monocle3轨迹分析
##创建CDS对象并预处理数据
scRNASeJeAP=readRDS('../cellchat/scRNASeJeAPharmony.rda')

# SeuratObject$Cluster<- paste0('cluster',SeuratObject$seurat_clusters)
scRNASeJeAP@meta.data$celltype2<- paste0(scRNASeJeAP@meta.data$cellType,scRNASeJeAP@meta.data$seurat_clusters)

SeuratObject2=scRNASeJeAP

p1=DimPlot(SeuratObject2,group.by = "orig.ident")
# 查看聚类信息
p2=DimPlot(SeuratObject2,group.by = "seurat_clusters",label=T)
p3=DimPlot(SeuratObject2,group.by = "celltype2",label=T)
combine<-CombinePlots(list(p1,p2,p3))
ggsave(filename = "p_DimPlot_sample_seurat_clusters.pdf", plot = combine, height = 6, width = 7)


SeuratObject=SeuratObject2
data <- GetAssayData(SeuratObject, assay = 'RNA', slot = 'counts')
cell_metadata <- SeuratObject@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 20)
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
ggsave(filename = "cds_int.umap.pdf", plot = combine, height = 4, width = 7)
## Monocle3聚类分区
cds <- cluster_cells(cds,resolution = 0.004)
## 识别轨迹
cds <- learn_graph(cds)

p1 <- plot_cells(cds, color_cells_by='celltype2',show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
# p = wrap_plots(p1, p2)


p3 = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = T, color_cells_by='celltype2', 
               label_branch_points = T)
# p4 = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = F, color_cells_by='celltype2', 
#                label_branch_points = F)
# ggsave(filename = "cds_int.umap-p4.pdf", plot = p4, height = 4, width = 7)
##细胞按拟时排序
##############################################################
# cds <- order_cells(cds) 存在bug，使用辅助线选择root细胞
options(browser="firefox")# 
# 通过在终端中执行以下命令来设置 DISPLAY 环境变量
# export DISPLAY=:0.0

# utils::browseURL("http://google.com/", browser = '/usr/bin/firefox')交互操作。
# utils::browseURL("http://baidu.com/", browser = '/usr/bin/firefox')

cds <- order_cells(cds)
p4=                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         (cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE,  label_branch_points = FALSE)
combine<-CombinePlots(list(p1,p2,p3,p4))
ggsave(filename = "cds_int.umap2.pdf", plot = combine, height = 6, width = 7)



Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
#挑选top10画图展示
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#基因表达趋势图
pdf(file="cds_gene_track_top10.pdf",width=6.5,height=6)
plot_genes_in_pseudotime(cds[Track_genes_sig,], 
                         min_expr=0.5, ncol = 2)
dev.off()
#FeaturePlot图
pdf(file="cds_gene_featureplot_top10.pdf",width=6.5,height=6)
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
dev.off()
##寻找共表达模块
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)$predicted.id)
cell_group[,2]=cds@colData[,'celltype2']
agg_mat <- aggregate_gene_expression(cds,gene_group_df=gene_module, cell_group_df=cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

pdf(file="cds_gene_top10hm.pdf",width=4,height=12)
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
dev.off()



# 
cds <- readRDS('./m3_cds.rds')
class(cds)
saveRDS(cds,file='m3_cds_order.rds')
test_genes <- c("Nucb2","Igf1r")
test_lineage_cds <- cds[rowData(cds)$gene_short_name %in% test_genes,]
test_lineage_cds <- order_cells(test_lineage_cds)

pdf(file="cds_gene_pseu.pdf",width=4,height=4)
plot_genes_in_pseudotime(test_lineage_cds,
                         color_cells_by="celltype2",
                         min_expr=0.5)
dev.off()
pdf(file="cds_gene_celltype.pdf",width=4,height=4)
plot_genes_by_group(cds,
                    test_genes,
                    group_cells_by="celltype2",
                    ordering_type="maximal_on_diag",
                    max.size=3)
dev.off()



cds <- readRDS('m3_cds_order.rds')
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
#挑选top10画图展示
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#基因表达趋势图
pdf(file="cds_gene_track_top10.pdf",width=6.5,height=6)
plot_genes_in_pseudotime(cds[Track_genes_sig,], 
                         min_expr=0.5, ncol = 2)
dev.off()



pdf(file="cds_gene_track_NUCB2IGF1R.pdf",width=6.5,height=6)
plot_genes_in_pseudotime(cds[c('Nucb2','Igf1r'),], 
                         min_expr=0.5, ncol = 2)
dev.off()
