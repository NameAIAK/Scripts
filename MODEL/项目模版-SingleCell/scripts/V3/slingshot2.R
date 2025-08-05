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
scale.data <- VCM@assays$Spatial@scale.data
scale.gene <- rownames(scale.data)

counts <- VCM@assays$Spatial@counts


# 将表达矩阵转换为SingleCellExperiment对象
# 输入需要counts矩阵，否则影响下游分析
sim <- SingleCellExperiment(assays = List(counts = counts)) 

# umap reduction
umap = VCM@reductions$umap@cell.embeddings
colnames(umap) = c('UMAP-1', 'UMAP-2')
# 将降维的结果添加到SingleCellExperiment对象中
reducedDims(sim) = SimpleList(UMAP = umap)

# metadata
meta = VCM@meta.data
# colData(sim)相当于meta.data，但他不是data.frame格式
# 所以需要分步赋予
colData(sim)$sampleId = meta$orig.ident
colData(sim)$group = meta$group
colData(sim)$subtype = meta$subtype

rd = umap
plot(rd, col = rgb(0,0,0,.5), pch=16, asp = 1)
# 一行命令就可生成，这里计算是比较快的
sim <- slingshot(sim, 
                 clusterLabels = 'subtype',  # 选择colData中细胞注释的列名
                 reducedDim = 'UMAP',  
                 start.clus= NULL,  # 选择起点
                 end.clus = NULL     # 这里我们不指定终点
)     
colnames(colData(sim))
plot(reducedDims(sim)$UMAP, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(9,"Set1"))
legend("bottomleft",
       legend = paste0("lineage",1),
       col = unique(brewer.pal(6,"Set1")),
       inset=0.8,
       pch = 16)


summary(sim$slingPseudotime_1)


# 我们做一个绚丽的渐变色彩虹色
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) # 我们把这些颜色变成100个梯度，模拟渐变色
plotcol <- colors[cut(-sim$slingPseudotime_1, breaks=100)] # 这里我们用cut函数把 lineage1 分割成100个区间，同一区间的细胞给与同一个颜色
plotcol[is.na(plotcol)] <- "lightgrey" # 不属于 lineage3 的点为NA，我们把他们变成灰色
plotcol


plot(reducedDims(sim)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(9,"Set1"))
legend("bottomleft",
       legend = paste0("lineage",1),
       col = unique(brewer.pal(6,"Set1")),
       inset=0.8,
       pch = 16)

VCM$pseu_time<- as.numeric(sim$slingPseudotime_1)
SpatialFeaturePlot(VCM,features = 'pseu_time',pt.size.factor = 0.8)+
  scale_colour_gradientn(colours = plotcol)
FeaturePlot(ne1,features ='pseu_time')+
  scale_colour_gradientn(colours = plotcol)




















