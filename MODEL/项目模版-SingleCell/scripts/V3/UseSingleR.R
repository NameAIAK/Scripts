library(SingleR)
library(celldex)
library(Seurat)
# 使用SingleR的最简单方法是使用内置参考对细胞进行注释。singleR自带的7个参考数据集，其中5个是人类数据，2个是小鼠的数据：
# BlueprintEncodeData Blueprint (Martens and Stunnenberg 2013) and Encode (The ENCODE Project Consortium 2012) （人）
# DatabaseImmuneCellExpressionData The Database for Immune Cell Expression(/eQTLs/Epigenomics)(Schmiedel et al. 2018)（人）
# HumanPrimaryCellAtlasData the Human Primary Cell Atlas (Mabbott et al. 2013)（人）
# MonacoImmuneData, Monaco Immune Cell Data - GSE107011 (Monaco et al. 2019)（人）
# NovershternHematopoieticData Novershtern Hematopoietic Cell Data - GSE24759（人）
# ImmGenData the murine ImmGen (Heng et al. 2008) （鼠）
# MouseRNAseqData a collection of mouse data sets downloaded from GEO (Benayoun et al. 2019).鼠）

##下载注释数据库
# BlueprintEncodeData
# DatabaseImmuneCellExpressionData
# HumanPrimaryCellAtlasData
# ImmGenData
# MonacoImmuneData
# MouseRNAseqData
# NovershternHematopoieticData

hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
save(hpca.se,file='HumanPrimaryCellAtlas.RData')
#直接load下载好的数据库
load("HumanPrimaryCellAtlas.RData")
# load("BlueprintEncode.RData")

# 加载单细胞数据
scRNA=readRDS("../ana/scRNA.celltype.rds")

# SingleR鉴定细胞类型
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = hpca.se, labels = hpca.se$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
# 使用多个数据库注释
# ref = list(BP=bpe.se, HPCA=hpca.se), labels = list(bpe.se$label.main, hpca.se$label.main)  

cellType = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(cellType,"celltype_singleR.csv",row.names = F)
scRNA@meta.data$cellType = "NA"
for(i in 1:nrow(cellType)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == cellType$ClusterID[i]),'cellType'] <- cellType$celltype[i]}

# 鉴定结果展示
library(patchwork)
p1 = DimPlot(scRNA, group.by="cellType", label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNA, group.by="cellType", label=T, label.size=5, reduction='umap')
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')

ggsave("tSNE_celltype.pdf", p1, width=7 ,height=6)
ggsave("UMAP_celltype.pdf", p2, width=7 ,height=6)
ggsave("celltype.pdf", p3, width=10 ,height=5)


# 总结：结果不建议使用，不太行（抠鼻）