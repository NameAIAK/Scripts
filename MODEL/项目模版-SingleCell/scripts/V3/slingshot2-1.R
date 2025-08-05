# 首先你需要把这些包安装好
suppressPackageStartupMessages({
  library(slingshot)
  library(SingleCellExperiment)
  library(qs)
  library(tidyverse)
  library(RColorBrewer)
})
library(Seurat)
# 在这里为了减少计算压力，我们把高变基因提取出来
# 你可以跟我一样从 scale.data 中提取，也可以直接在 seurat 对象中找出来
scRNASeJeAP=readRDS('../cellchat/scRNASeJeAPharmony.rda')

# SeuratObject$Cluster<- paste0('cluster',SeuratObject$seurat_clusters)
scRNASeJeAP@meta.data$celltype2<- paste0(scRNASeJeAP@meta.data$cellType,scRNASeJeAP@meta.data$seurat_clusters)
scale.data <- scRNASeJeAP@assays$RNA@scale.data
scale.gene <- rownames(scale.data)

counts <- scRNASeJeAP@assays$RNA@counts


# 将表达矩阵转换为SingleCellExperiment对象
# 输入需要counts矩阵，否则影响下游分析
sim <- SingleCellExperiment(assays = List(counts = counts)) 

# umap reduction
umap = scRNASeJeAP@reductions$umap@cell.embeddings
colnames(umap) = c('UMAP-1', 'UMAP-2')
# 将降维的结果添加到SingleCellExperiment对象中
reducedDims(sim) = SimpleList(UMAP = umap)

# metadata
meta = scRNASeJeAP@meta.data
# colData(sim)相当于meta.data，但他不是data.frame格式
# 所以需要分步赋予
colData(sim)$sampleId = meta$orig.ident
# colData(sim)$group = meta$group
colData(sim)$subtype = meta$celltype2

rd = umap
pdf(file="slingshot_rd.pdf",width=4,height=3.5)
plot(rd, col = rgb(0,0,0,.5), pch=16, asp = 1)
dev.off()
# 一行命令就可生成，这里计算是比较快的
sim <- slingshot(sim, 
                 clusterLabels = 'subtype',  # 选择colData中细胞注释的列名
                 reducedDim = 'UMAP',  
                 start.clus= NULL,  # 选择起点
                 end.clus = NULL     # 这里我们不指定终点
)     
colnames(colData(sim))

pdf(file="slingshot_umapline_single.pdf",width=4,height=3.5)
plot(reducedDims(sim)$UMAP, pch=16, asp = 1)
lines(SlingshotDataSet(sim[[1]]), lwd=2, col=brewer.pal(9,"Set1"))
legend("bottomleft",
       legend = paste0("lineage",1),
       col = unique(brewer.pal(6,"Set1")),
       inset=0.8,
       pch = 16)
dev.off()

summary(sim$slingPseudotime_1)


# 我们做一个绚丽的渐变色彩虹色
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) # 我们把这些颜色变成100个梯度，模拟渐变色
plotcol <- colors[cut(-sim$slingPseudotime_1, breaks=100)] # 这里我们用cut函数把 lineage1 分割成100个区间，同一区间的细胞给与同一个颜色
plotcol[is.na(plotcol)] <- "lightgrey" # 不属于 lineage1 的点为NA，我们把他们变成灰色
# plotcol

pdf(file="slingshot_umapline.pdf",width=7,height=5.5)
plot(reducedDims(sim)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(9,"Set1"))
legend("topleft",
       legend = paste0("lineage",1:3),
       col = unique(brewer.pal(6,"Set1")),
       inset=0.8,
       pch = 16)
dev.off()

pdf(file="slingshot_umapline_single.pdf",width=7,height=5.5)
plot(reducedDims(sim)$UMAP,  col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim), linInd = 1, lwd=2, col=brewer.pal(9,"Set1"))
legend("topleft",
       legend = paste0("lineage",1),
       col = unique(brewer.pal(6,"Set1")),
       inset=0.8,
       pch = 16)
dev.off()

#-----------------------------------step5 tradeSeq 下游分析
library(tradeSeq)
# Fit negative binomial model
counts <- sim@assays@data$counts
crv <- SlingshotDataSet(sim)
#拟合负二项式模型 需要决定结的数量
set.seed(111)
icMat <- evaluateK(counts = counts, 
                   sds = crv, 
                   k = 3:10,    # no more than 12
                   nGenes = 500,
                   verbose = T)
 
set.seed(111)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
# fit negative binomial GAM
# 2k cells ~13 min
# system.time()这个函数可以计算运行时间
 
sce <- fitGAM(counts = counts, 
                      pseudotime = pseudotime, 
                      cellWeights = cellWeights,
                      nknots = 6, 
                      verbose = FALSE)
 
#探索基因表达与拟时序的相关性
assoRes <- associationTest(sce)
head(assoRes)
#寻找与起止点相关性最高的基因
startRes <- startVsEndTest(sce)
head(startRes)
# 按相关性进行排序
oStart <- order(startRes$waldStat, decreasing = TRUE)
# 挑选相关性最强的基因，并可视化
sigGeneStart <- names(sce)[oStart[1]]
pdf(file="slingshot_plotsmoothers.pdf",width=4,height=4)
plotSmoothers(sce, counts, gene = sigGeneStart)
dev.off()
#--------------美化基因表达量变化图
# 取celltpye和配色信息
coldata <- data.frame(celltype = sim@colData$subtype,
                      plotcol = plotcol)
rownames(coldata) = colnames(sim)
 
# 把sce中的3000个细胞对应信息取出
filter_coldata <- coldata[colnames(sce),]
 
# 添加拟时序信息
filter_coldata$Pseudotime = sce$crv$pseudotime.Lineage1
 
# top6 genes
top6 <- names(sce)[oStart[1:6]]
top6_exp = sce@assays@data$counts[top6,] 
top6_exp = log2(top6_exp + 1) %>% t()
 
# 获得最终数据
plt_data = cbind(filter_coldata, top6_exp)
colnames(plt_data)
 
 
#画图
#画图
library(ggplot2)
library(ggsci)
library(ggpubr)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
mycolors = getPalette(length(unique(top6))+2)
 
# 为了拼图美观，我们把legend隐藏掉
plt_list = list()
for (gene in top6) {
        print(gene)
        p = ggscatter(data = plt_data,
                      x = 'Pseudotime',
                      y = gene,
                      color = 'celltype',
                      size = 0.6)+
                geom_smooth(se = F, color = 'orange')+
                theme_bw()+
                scale_color_manual(values = mycolors)+
                theme(legend.position = 'right')
        plt_list[[gene]] = p
}
 
library(patchwork)
library(ggpubr)
library(Rmisc)
grid_plot_up <- plot_grid(plotlist =plt_list , ncol = 3)
# grid_plot_up
ggsave(
  plot = grid_plot_up,
  "test.pdf",
  height = 12,
  width = 16
)

# top6 <- names(sce)[oStart[1:6]]
gene_exp = sce@assays@data$counts[c('Nucb2','Igf1r'),] 
gene_exp = log2(gene_exp + 1) %>% t()
 
# 获得最终数据
plt_data2 = cbind(filter_coldata, gene_exp)
# colnames(plt_data)
p = ggscatter(data = plt_data2,
              x = 'Pseudotime',
              y = 'Nucb2',
              color = 'celltype',
              size = 0.6)+
       geom_smooth(se = F, color = 'orange')+
       theme_bw()+
       scale_color_manual(values = mycolors)+
       theme(legend.position = 'right')
ggsave(plot=p,'Nucb2.pdf')

p = ggscatter(data = plt_data2,
              x = 'Pseudotime',
              y = 'Igf1r',
              color = 'celltype',size = 0.6)+
       geom_smooth(se = F, color = 'orange')+
       theme_bw()+
       scale_color_manual(values = mycolors)+
       theme(legend.position = 'right')
ggsave(plot=p,'Igf1r.pdf')

sce_slingshot1=sim
mouse_data=scRNASeJeAP
sce_slingshot1$sling_pseudotime = sce_slingshot1[[paste0("slingPseudotime_1")]]
mouse_data$sling_pseudotime = sce_slingshot1$sling_pseudotime

#仅分析在轨迹拟时中的细胞，去除NA
sce_slingshot1_l1 = sce_slingshot1[,!is.na(sce_slingshot1$sling_pseudotime)]
seur = mouse_data[,!is.na(mouse_data$sling_pseudotime)]

#和前面一样这一步会比较慢，这个数据跑了大概1个小时。
sce_slingshot1_l1 <- fitGAM(counts(sce_slingshot1_l1), 
                            cellWeights = rep(1, ncol(sce_slingshot1_l1)), 
                            pseudotime = sce_slingshot1_l1$sling_pseudotime)

ATres <- associationTest(sce_slingshot1_l1)
association_test_tab = as_tibble(cbind(gene = rownames(ATres), ATres))
slingshot_for_plotMatrix <- function(seurat_obj,
                                     n_bins,#拟时需要分割的区间，将相似的拟时区间合并，这类似于我们monocle3中的方式
                                     min_exp)
                                     {
  seurat_meta = seurat_obj@meta.data
  seurat_meta = as_tibble(cbind(cell.id = as.character(rownames(seurat_meta)), seurat_meta))
  seurat_meta = seurat_meta[order(seurat_meta$sling_pseudotime),]
  
  pl_cells = as.character(seurat_meta$cell.id)
  
  #提取表达矩阵,并将cell id的exp排序与前面排序好的cell id一致
  exp = seurat_obj@assays$RNA@data
  exp = exp[,colnames(exp) %in% pl_cells]
  expr_mat = exp[,order(match(colnames(exp), pl_cells))]
  
  expr_mat = as.matrix(expr_mat[rownames(expr_mat) %in% association_test_tab$gene,])
  
  clust_expr_mat = matrix(nrow = nrow(expr_mat), 
                          ncol = n_bins, dimnames = list(rownames(expr_mat), 1:n_bins))
  
  max_pseudotime = max(seurat_meta$sling_pseudotime)
  pseudotime_bin_size = max_pseudotime/n_bins
  
  pseudotime_cluster_stat = NULL
  seurat_obj$pseudotime_bin = NA_integer_
  
  for (i in 1 : n_bins){
    
    bin_cells = seurat_meta$cell.id[(seurat_meta$sling_pseudotime > (i-1)*pseudotime_bin_size & 
                                       seurat_meta$sling_pseudotime <= i*pseudotime_bin_size)]

    
    seurat_obj$pseudotime_bin[colnames(seurat_obj) %in% bin_cells] = i

    #计算基因平均表达量
    if (length(bin_cells)>10){
      m2 = expr_mat[,colnames(expr_mat) %in% bin_cells]
      clust_expr_mat[,i] = apply(m2, 1, mean, na.rm = TRUE)
    }

  }
  
  #数据缩放一下，为了更好的展现热图，并删除低表达基因
  mm1 = clust_expr_mat - apply(clust_expr_mat, 1, mean, na.rm = TRUE)
  mm2 = mm1[apply(abs(mm1),1, max, na.rm = TRUE)>min_exp,]
  
  return(mm2)
  
 }
mm = slingshot_for_plotMatrix(seurat_obj = seur, n_bins = 20, min_exp = 0.2)

max_range = max(range(is.finite(mm)))
lim = c(-max_range, max_range)

library(pheatmap)
heatmap1 = pheatmap(mm, show_rownames=F, cluster_rows = TRUE,
                    cluster_cols = FALSE, show_colnames = FALSE, 
                    clustering_distance_rows = "euclidean",
                    clustering_method = "ward.D2",
                    treeheight_row = 50,
                    cutree_rows = 5, 
                    color = colorRampPalette(rev(brewer.pal(9, "PRGn")))(250),
                    breaks = seq(lim[1]/4, lim[2]/4, length.out = 251),
                    border_color = NA)

# ##################
annotation_row <- data.frame(Cluster=factor(cutree(heatmap1$tree_row, 5)))
row.names(annotation_row) <- rownames(mm)


rowcolor <- c("#85B22E","#E29827","#922927",'#57C3F3','#A758E5')
names(rowcolor) <- c("1","2","3","4","5") #类型颜色

ann_colors <- list(pseudotime=viridis(100),
                  Cluster=rowcolor) #颜色设置
                    
pheatmap(mm,
        useRaster = T,
        cluster_cols=F,
        cluster_rows=T,
        show_rownames=F,
        show_colnames=F,
        clustering_method = "ward.D2",
        clustering_distance_rows = "euclidean",
        cutree_rows=5,
        border_color = NA,
        filename=NA,
        color = colorRampPalette(rev(brewer.pal(9, "PRGn")))(100),
        breaks = seq(lim[1], lim[2], length.out = 100),
        annotation_colors=ann_colors,
        annotation_row = annotation_row
)






























# 空间单细胞
# scRNASeJeAP$pseu_time<- as.numeric(sim$slingPseudotime_1)

# p1=SpatialFeaturePlot(scRNASeJeAP,features = 'pseu_time',pt.size.factor = 0.8)+
#   scale_colour_gradientn(colours = plotcol)
# p2=FeaturePlot(ne1,features ='pseu_time')+
#   scale_colour_gradientn(colours = plotcol)
# combine=CombinePlots(list(p1,p2))
# ggsave(filename = "slingshot_pseu.pdf", plot = combine, height = 6, width = 7)




















