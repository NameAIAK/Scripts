###安装
install.packages("BiocManager")
BiocManager::install("monocle")
install.packages("RColorBrewer") #画图配色板
# BiocManager::install("monocle")
# BiocManager::install("HSMMSingleCell")
# install.packages("F:\\R包\\HSMMSingleCell_1.24.0.tar.gz", repos = NULL, type = "source")
library(monocle)
library(RColorBrewer)
library(Seurat)
library(tidyverse)
library(patchwork)

# remove.packages("igraph")
# require(devtools)
# install_version("igraph", version = "2.0.3",repos = "http://cran.us.r-project.org")
# install.packages('../igraph_2.0.3.tar.gz',repos=NULL,type='source')

###数据准备
##mel.R
scRNA <- readRDS("../scRNA_seurat_cellType.RDS")
scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))
# scRNASeJeAPmyoFBECFB=subset(x = scRNASeJe, subset = (cellType == "AP"|cellType == "myoFB"|cellType == "EC"|cellType == "FB"))
# SeuratObject=scRNASeJeAPmyoFBECFB
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

# ###创建对象
# # 创建monocle对象
# cd <- newCellDataSet(as.matrix(expr_matrix), 
#                      phenoData = pd, 
#                      featureData = fd, 
#                      lowerDetectionLimit = 0.1, 
#                      expressionFamily = tobit(Lower = 0.1))
# # expressionFamily: 数据为TPM/FPKM时设置为tobit(Lower = 0.1)，数据为count时设置为negbinomial.size())

# # 将FPKM/TPM数据转换为UMI数据（read count）
# rpc_matrix <- relative2abs(cd)

# 重新创建monocle对象
# expr_matrix
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
# 设置颜色
color1 <- c(brewer.pal(6, "Set1"))
getPalette <- colorRampPalette(brewer.pal(6, "Set1"))
# 按照表型进行映射
pdf(file="cd_cellType.trajectory.pdf",width=6.5,height=6)
plot_cell_trajectory(cd, show_cell_names = F, color_by = "celltype2") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"celltype2"]))))
dev.off()

# 按照阶段进行映射（手动设置颜色）
pdf(file="cd_cluster.trajectory.pdf",width=6.5,height=6)
plot_cell_trajectory(cd, show_cell_names = F, color_by = "Cluster") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"Cluster"]))))
dev.off()

# 按照计算的Pseudotime进行映射
pdf(file="cd_Pseudotime.trajectory.pdf",width=6.5,height=6)
plot(plot_cell_trajectory(cd, show_cell_names = F, color_by = "Pseudotime") + scale_color_viridis_c())
dev.off()

# 按照基因表达进行映射
# Igf1r,Nebl
cd$exprs<- SeuratObject@assays$RNA@counts['Nucb2',]
pdf(file="cd_gene_Nucb2.trajectory.pdf",width=6.5,height=6)
plot(plot_cell_trajectory(cd, show_cell_names = F, color_by = "exprs") + scale_color_viridis_c())
dev.off()

sig_genes <- subset(diff_test_res, qval < 1e-5)
sig_genes<- sig_genes[order(sig_genes$qval),]
cg<- as.character(head(sig_genes$gene_short_name))

# 画差异基因的表达图
pdf(file="cd_gene_expression.trajectory.pdf",width=6.5,height=6)
plot_genes_jitter(cd[cg,],
                  grouping = 'cellType',
                  color_by = 'cellType',
                  nrow = 3,
                  ncol = NULL)
dev.off()


pData(cd)$sample=SeuratObject$orig.ident
nrows_state=1
pdf(file="cd_celltype2.trajectory.pdf",width=2,height=3.5)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"cellType"]))))
plot_cell_trajectory(cd, color_by = "celltype2", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~sample, nrow = nrows_state, scales = "free")+ scale_color_manual(values=getPalette(length(unique(sample_sheet[,"celltype2"]))))+ggtitle("cellType") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()

pdf(file="cd_cellType.trajectory.pdf",width=6.5,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"cellType"]))))
plot_cell_trajectory(cd, color_by = "cellType", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~sample, nrow = nrows_state, scales = "free")+ scale_color_manual(values=getPalette(length(unique(sample_sheet[,"cellType"]))))+ggtitle("cellType") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()

pdf(file="cd_sample.trajectory.pdf",width=6.5,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"cellType"]))))
plot_cell_trajectory(cd, color_by = "sample", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~cellType, nrow = nrows_state, scales = "free")+ggtitle("sample") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()

pdf(file="cd_sample2.trajectory.pdf",width=10,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "cellType") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"cellType"]))))
plot_cell_trajectory(cd, color_by = "sample", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~celltype2, nrow = nrows_state, scales = "free")+ggtitle("sample") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()

# 按照阶段进行映射（手动设置颜色）
pdf(file="cd_cluster.trajectory.pdf",width=6.5,height=6)
# plot_cell_trajectory(cd, show_cell_names = F, color_by = "Cluster") + scale_color_manual(values = getPalette(length(unique(sample_sheet[,"Cluster"]))))
plot_cell_trajectory(cd, color_by = "Cluster", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~sample, nrow = nrows_state, scales = "free")+ scale_color_manual(values=getPalette(length(unique(sample_sheet[,"Cluster"]))))+ggtitle("Cluster") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()

# 按照计算的Pseudotime进行映射
pdf(file="cd_Pseudotime.trajectory.pdf",width=6.5,height=6)
# plot(plot_cell_trajectory(cd, show_cell_names = F, color_by = "Pseudotime") + scale_color_viridis_c())
plot_cell_trajectory(cd, color_by = "Pseudotime", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~sample, nrow = nrows_state, scales = "free")+ggtitle("Pseudotime") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()

cd$exprs<- SeuratObject@assays$RNA@counts['Nucb2',]
pdf(file="cd_gene_Nucb2.trajectory.pdf",width=6.5,height=6)
# plot(plot_cell_trajectory(cd, show_cell_names = F, color_by = "exprs") + scale_color_viridis_c())
plot(plot_cell_trajectory(cd, color_by = "exprs", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~sample, nrow = nrows_state, scales = "free")+scale_color_viridis_c())
dev.off()

cd$exprs2<- SeuratObject@assays$RNA@counts['Igf1r',]
pdf(file="cd_gene_Igf1r.trajectory.pdf",width=6.5,height=6)
# plot(plot_cell_trajectory(cd, show_cell_names = F, color_by = "exprs") + scale_color_viridis_c())
plot(plot_cell_trajectory(cd, color_by = "exprs2", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~sample, nrow = nrows_state, scales = "free")+scale_color_viridis_c())
dev.off()

#########plan2#########
###安装
install.packages("BiocManager")
BiocManager::install("monocle")
install.packages("RColorBrewer") #画图配色板
library(monocle)
library(RColorBrewer)
library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)

###数据准备
##mel.R
scRNA <- readRDS("../scRNA_seurat_cellType.RDS")
scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))
scRNASeJeAPmyoFB=subset(x = scRNASeJe, subset = (cellType == "AP"))
SeuratObject=scRNASeJeAPmyoFB
scRNASeJeAP=readRDS('../cellchat/scRNASeJeAPharmony.rda')

scRNASeJeAP@meta.data$cellType2<- paste0(scRNASeJeAP@meta.data$cellType,scRNASeJeAP@meta.data$seurat_clusters)
 
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

# SeuratObject=subset(x = SeuratObject, subset = (celltype2 == "AP3"|celltype2 == "AP6"))

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

# 创建monocle对象
# expr_matrix
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
saveRDS(cd,file='m2_cd_order.rds')
# saveRDS(cd,file='m2_cd_order_sub.rds')
# 排序好的细胞可以进行可视化，可标注细胞的各注释信息)
# 查看细胞注释信息
head(cd@phenoData@data)
# 设置颜色
color1 <- c(brewer.pal(6, "Set1"))
getPalette <- colorRampPalette(brewer.pal(6, "Set1"))

cd <- readRDS('m2_cd_order.rds')
to_be_tested <- row.names(subset(fData(cd),
                                 gene_short_name %in% c("Nucb2","Igf1r")))
cds_subset <- cd[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
# diff_test_res[,c("gene_short_name", "pval", "qval")]
pdf(file="m2_genes_in_pseudotime.pdf",width=4,height=3)
plot_genes_in_pseudotime(cds_subset, color_by = "celltype2")
dev.off()

diff_test_res <- differentialGeneTest(cd[to_be_tested,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(cd[sig_gene_names,],
                        num_clusters = 6,
                        cores = 1,
                        show_rownames = T)
pdf(file="m2ap07celltype2.pdf",height=6,width=6)
plot_cell_trajectory(cd, color_by = "celltype2")
dev.off()
# 图如下
BEAM_res <- BEAM(cd, branch_point = 3, cores = 10, progenitor_method = "duplicate")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
pdf(file="m2_genes_branched_heatmap.pdf",width=4,height=3)
plot_genes_branched_heatmap(cd[c(row.names(subset(BEAM_res[1:20,],qval < 1e-4)),'Nucb2','Igf1r'),],
                            branch_point = 3,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()


cd_sub=readRDS('./m2_cd_order_sub.rds')
to_be_tested2 <- row.names(subset(fData(cd_sub),
                                 gene_short_name %in% c("Nucb2","Igf1r")))
cds_subset <- cd_sub[to_be_tested2,]
diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
# diff_test_res[,c("gene_short_name", "pval", "qval")]
pdf(file="m2sub_genes_in_pseudotime.pdf",width=4,height=3)
plot_genes_in_pseudotime(cds_subset, color_by = "celltype2")
dev.off()

