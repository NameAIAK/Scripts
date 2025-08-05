library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(Seurat)
library(SeuratData)


cells_rankings <- AUCell_buildRankings(SObj@assays$RNA@data, nCores=1, plotStats=TRUE) 
cells_rankings
#download:https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C5
c5 <- read.gmt("C:/Users/LXF/Desktop/c5.go.cc.v7.4.symbols.gmt")
head(c5$term)
head(c5$gene)
geneSets <- lapply(unique(c5$term), function(x){print(x);c5$gene[c5$term == x]})
names(geneSets) <- unique(c5$term)
geneSets$GOCC_SAGA_COMPLEX

?AUCell_calcAUC
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
cells_AUC
length(rownames(cells_AUC@assays@data$AUC))

grep("REG",rownames(cells_AUC@assays@data$AUC),value = T)

##set gene set of interest here for plotting
# aucs <- as.numeric(getAUC(cells_AUC["GOCC_NUCLEOTIDE_EXCISION_REPAIR_COMPLEX",]))
aucs <- getAUC(cells_AUC)
scRNASeJe$AUC <- aucs
#用mel.R的结果
df<- data.frame(scRNASeJe@meta.data, scRNASeJe@reductions$umap@cell.embeddings)
head(df)

class_avg <- df %>%
  group_by(seurat_clusters) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

pdf('gsva.pdf',,height = 6, width = 6)
ggplot(df, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="A") +
  ggrepel::geom_label_repel(aes(label = seurat_clusters),
                            data = class_avg,
                            size = 6,
                            label.size = 0,
                            segment.color = NA
  )+   theme(legend.position = "none",title='aa') + theme_bw()+labs(title='mitosis')
dev.off()
?theme

####Se 通路活性可视化
#download:https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C5
fatty <- read.gmt("fatty.gmt")
head(fatty$term)
head(fatty$gene)
geneSets <- lapply(unique(fatty$term), function(x){print(x);fatty$gene[fatty$term == x]})
names(geneSets) <- unique(fatty$term)

cells_rankings <- AUCell_buildRankings(scRNASe@assays$RNA@data, plotStats=TRUE) 
cells_rankings
# geneSets$GOCC_SAGA_COMPLEX

# ?AUCell_calcAUC
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
# cells_AUC
# length(rownames(cells_AUC@assays@data$AUC))

# grep("FATTY",rownames(cells_AUC@assays@data$AUC),value = T)

# 循环作各通路dimplot
p_list <- list() 
i=1
for (pathway in rownames(cells_AUC@assays@data$AUC)){
    aucs <- as.numeric(getAUC(cells_AUC[pathway,]))
    # Idents(scRNASe)='cluster'
    scRNASe$AUC <- aucs
    #用mel.R的结果
    df<- data.frame(scRNASe@meta.data, scRNASe@reductions$umap@cell.embeddings)
    head(df)

    class_avg <- df %>%
    group_by(cellType) %>%
    summarise(
        UMAP_1 = median(UMAP_1),
        UMAP_2 = median(UMAP_2)
    )

    
    p=ggplot(df, aes(UMAP_1, UMAP_2))  +
    geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="A") +
    ggrepel::geom_label_repel(aes(label = cellType),
                                data = class_avg,
                                size = 6,
                                label.size = 0,
                                segment.color = NA
    )+   theme(legend.position = "none",title='aa') + theme_bw()+labs(title=pathway)
    
    p_list[[i]] <- p
    i=i+1
}

grid_plot_up <- plot_grid(plotlist =p_list , nrow = 4,ncol = 5)

ggsave(
  plot = grid_plot_up,
  "Se.pdf",
  height = 40,
  width = 40
)




##################
fatty <- read.gmt("fat.gmt")
head(fatty$term)
head(fatty$gene)
geneSets <- lapply(unique(fatty$term), function(x){print(x);fatty$gene[fatty$term == x]})
names(geneSets) <- unique(fatty$term)
features=rownames(unique(fatty$term))
features=gsub( "_", "-",features)

############scRNASe
cells_rankings <- AUCell_buildRankings(scRNASe@assays$RNA@data, plotStats=TRUE) 
cells_rankings

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- getAUC(cells_AUC)
# aucs=t(aucs)
scRNASe_AUC <- CreateAssayObject(counts = aucs)
scRNASe@assays$AUC <- scRNASe_AUC
scRNASe@assays$AUC@key <- "rna_"
DefaultAssay(scRNASe) <- "AUC"

# 画Se通路
Idents(scRNASe)='cellType'
pdf(file="Se-Activepathway.pdf",height = 6, width = 18)
DotPlot(scRNASe,group.by = 'cellType', features = unique(features)) + coord_flip()+RotatedAxis()
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()

############scRNAJe
cells_rankings <- AUCell_buildRankings(scRNAJe@assays$RNA@data, plotStats=TRUE) 
cells_rankings

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- getAUC(cells_AUC)
# aucs=t(aucs)
scRNAJe_AUC <- CreateAssayObject(counts = aucs)
scRNAJe@assays$AUC <- scRNAJe_AUC
scRNAJe@assays$AUC@key <- "rna_"
DefaultAssay(scRNAJe) <- "AUC"

# 画Je通路
Idents(scRNAJe)='cellType'
pdf(file="Je-Activepathway.pdf",height = 6, width = 18)
DotPlot(scRNAJe,group.by = 'cellType', features = unique(features)) + coord_flip()+RotatedAxis()
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()

############scRNASi
cells_rankings <- AUCell_buildRankings(scRNASi@assays$RNA@data, plotStats=TRUE) 
cells_rankings

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- getAUC(cells_AUC)
# aucs=t(aucs)
scRNASi_AUC <- CreateAssayObject(counts = aucs)
scRNASi@assays$AUC <- scRNASi_AUC
scRNASi@assays$AUC@key <- "rna_"
DefaultAssay(scRNASi) <- "AUC"

# 画Si通路
Idents(scRNASi)='cellType'
pdf(file="Si-Activepathway.pdf",height = 6, width = 18)
DotPlot(scRNASi,group.by = 'cellType', features = unique(features)) + coord_flip()+RotatedAxis()
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()


############scRNAJi
cells_rankings <- AUCell_buildRankings(scRNAJi@assays$RNA@data, plotStats=TRUE) 
cells_rankings

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- getAUC(cells_AUC)
# aucs=t(aucs)
scRNAJi_AUC <- CreateAssayObject(counts = aucs)
scRNAJi@assays$AUC <- scRNAJi_AUC
scRNAJi@assays$AUC@key <- "rna_"
DefaultAssay(scRNAJi) <- "AUC"

# 画Ji通路
Idents(scRNAJi)='cellType'
pdf(file="Ji-Activepathway.pdf",height = 6, width = 18)
DotPlot(scRNAJi,group.by = 'cellType', features = unique(features)) + coord_flip()+RotatedAxis()
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()


############scRNAST
cells_rankings <- AUCell_buildRankings(scRNAST@assays$RNA@data, plotStats=TRUE) 
cells_rankings

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- getAUC(cells_AUC)
# aucs=t(aucs)
scRNAST_AUC <- CreateAssayObject(counts = aucs)
scRNAST@assays$AUC <- scRNAST_AUC
scRNAST@assays$AUC@key <- "rna_"
DefaultAssay(scRNAST) <- "AUC"

# 画ST通路
Idents(scRNAST)='cellType'
pdf(file="ST-Activepathway.pdf",height = 6, width = 18)
DotPlot(scRNAST,group.by = 'cellType', features = unique(features)) + coord_flip()+RotatedAxis()
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()

############scRNAJT
cells_rankings <- AUCell_buildRankings(scRNAJT@assays$RNA@data, plotStats=TRUE) 
cells_rankings

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- getAUC(cells_AUC)
# aucs=t(aucs)
scRNAJT_AUC <- CreateAssayObject(counts = aucs)
scRNAJT@assays$AUC <- scRNAJT_AUC
scRNAJT@assays$AUC@key <- "rna_"
DefaultAssay(scRNAJT) <- "AUC"

# 画JT通路
Idents(scRNAJT)='cellType'
pdf(file="JT-Activepathway.pdf",height = 6, width = 18)
DotPlot(scRNAJT,group.by = 'cellType', features = unique(features)) + coord_flip()+RotatedAxis()
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()


############scRNA allsample
cells_rankings <- AUCell_buildRankings(scRNA@assays$RNA@data, plotStats=TRUE) 
cells_rankings

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- getAUC(cells_AUC)
# aucs=t(aucs)
scRNA_AUC <- CreateAssayObject(counts = aucs)
scRNA@assays$AUC <- scRNA_AUC
scRNA@assays$AUC@key <- "rna_"
DefaultAssay(scRNA) <- "AUC"
scRNA@meta.data$orig.ident=factor(scRNA@meta.data$orig.ident,levels =c('Se','Je','Si','Ji','ST','JT'))
# 画scRNA通路
Idents(scRNA)='orig.ident'
pdf(file="scRNA-Activepathway.pdf",height = 6, width = 12)
DotPlot(scRNA,group.by = 'orig.ident', features = unique(features)) + coord_flip()+RotatedAxis()
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()

saveRDS(scRNA, file = "scRNA-Activepathway.RDS")

############scRNA allsample AP
cells_rankings <- AUCell_buildRankings(scRNAAP@assays$RNA@data, plotStats=TRUE) 
cells_rankings

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- getAUC(cells_AUC)
# aucs=t(aucs)
scRNAAP_AUC <- CreateAssayObject(counts = aucs)
scRNAAP@assays$AUC <- scRNAAP_AUC
scRNAAP@assays$AUC@key <- "rna_"
DefaultAssay(scRNAAP) <- "AUC"
scRNAAP@meta.data$orig.ident=factor(scRNAAP@meta.data$orig.ident,levels =c('Se','Je','Si','Ji','ST','JT'))
# 画scRNA通路
Idents(scRNAAP)='orig.ident'
pdf(file="scRNAAP-Activepathway.pdf",height = 6, width = 12)
DotPlot(scRNAAP,group.by = 'orig.ident', features = unique(features)) + coord_flip()+RotatedAxis()
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()

saveRDS(scRNAAP, file = "scRNAAP-Activepathway.RDS")






scRNAAP=subset(x = scRNA, subset = (cellType == "AP"))
scRNASeJeAP=subset(x = scRNAAP, subset = (orig.ident == "Se"|orig.ident == "Je"))

##################
fatty <- read.gmt("fat.gmt")
head(fatty$term)
head(fatty$gene)
geneSets <- lapply(unique(fatty$term), function(x){print(x);fatty$gene[fatty$term == x]})
names(geneSets) <- unique(fatty$term)
features=rownames(unique(fatty$term))
features=gsub( "_", "-",features)

############scRNA allsample SeJe
scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))

cells_rankings <- AUCell_buildRankings(scRNASeJe@assays$RNA@data, plotStats=TRUE) 
cells_rankings

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- getAUC(cells_AUC)
# aucs=t(aucs)
scRNASeJe_AUC <- CreateAssayObject(counts = aucs)
scRNASeJe@assays$AUC <- scRNASeJe_AUC
scRNASeJe@assays$AUC@key <- "rna_"
DefaultAssay(scRNASeJe) <- "AUC"
scRNASeJe@meta.data$orig.ident=factor(scRNASeJe@meta.data$orig.ident,levels =c('Se','Je'))

# 画scRNA通路
Idents(scRNASeJe)='orig.ident'
pdf(file="scRNASeJe-Activepathway.pdf",height = 6, width = 12)
DotPlot(scRNASeJe,group.by = 'orig.ident', features = unique(features)) + coord_flip()+labs(x="Pathways")
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()

library(limma)
data=scRNASeJe@assays$AUC@data
dd_CM<- as.data.frame(as.matrix(data))
colnames(dd_CM)
table(Idents(scRNASeJe))
group2<-  factor(c(rep('Je', 10107),rep('Se', 8838)))#顺序靠后的 vs 顺序靠前的
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
write.table(lrt, 'SevsJe_1.txt', sep = '\t', col.names = NA, quote = FALSE)
features=c("BIOCARTA-VOBESITY-PATHWAY","WANG-CLASSIC-ADIPOGENIC-TARGETS-OF-PPARG","VERNOCHET-ADIPOGENESIS","GOBP-FAT-PAD-DEVELOPMENT","GOBP-NEGATIVE-REGULATION-OF-BROWN-FAT-CELL-DIFFERENTIATION")
Idents(scRNASeJe)='orig.ident'
pdf(file="scRNASeJe-ActivepathwayFC0.4.pdf",height = 3, width = 10)
DotPlot(scRNASeJe,group.by = 'orig.ident', features = unique(features)) + coord_flip()+labs(x="Pathways")
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()



############scRNA allsample SiJi
scRNASiJi=subset(x = scRNA, subset = (orig.ident == "Si"|orig.ident == "Ji"))
cells_rankings <- AUCell_buildRankings(scRNASiJi@assays$RNA@data, plotStats=TRUE) 
cells_rankings

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- getAUC(cells_AUC)
# aucs=t(aucs)
scRNASiJi_AUC <- CreateAssayObject(counts = aucs)
scRNASiJi@assays$AUC <- scRNASiJi_AUC
scRNASiJi@assays$AUC@key <- "rna_"
DefaultAssay(scRNASiJi) <- "AUC"
scRNASiJi@meta.data$orig.ident=factor(scRNASiJi@meta.data$orig.ident,levels =c('Si','Ji'))
# 画scRNA通路
Idents(scRNASiJi)='orig.ident'
pdf(file="scRNASiJi-Activepathway.pdf",height = 6, width = 12)
DotPlot(scRNASiJi,group.by = 'orig.ident', features = unique(features)) + coord_flip()+labs(x="Pathways")
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()


library(limma)
data=scRNASiJi@assays$AUC@data
dd_CM<- as.data.frame(as.matrix(data))
# colnames(dd_CM)
table(Idents(scRNASiJi))

# df <- df %>%
#   mutate(Group = case_when(
#     grepl("Si", Value) ~ "Si",
#     grepl("Ji", Value) ~ "Ji",
#     TRUE ~ "Other" # 处理不匹配任何上述条件的情况
#   ))
 
group2<-  factor(c(rep('Ji', 10025),rep('Si', 5800)))#顺序靠后的 vs 顺序靠前的
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
write.table(lrt, 'SivsJi_1.txt', sep = '\t', col.names = NA, quote = FALSE)


features=c("BIOCARTA-VOBESITY-PATHWAY","WANG-CLASSIC-ADIPOGENIC-TARGETS-OF-PPARG","VERNOCHET-ADIPOGENESIS","GOBP-NEGATIVE-REGULATION-OF-BROWN-FAT-CELL-DIFFERENTIATION")
Idents(scRNASiJi)='orig.ident'
pdf(file="scRNASiJi-ActivepathwayFC0.4.pdf",height = 3, width = 10)
DotPlot(scRNASiJi,group.by = 'orig.ident', features = unique(features)) + coord_flip()+labs(x="Pathways")
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()

############scRNA allsample STJT
scRNASTJT=subset(x = scRNA, subset = (orig.ident == "ST"|orig.ident == "JT"))
cells_rankings <- AUCell_buildRankings(scRNASTJT@assays$RNA@data, plotStats=TRUE) 
cells_rankings

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- getAUC(cells_AUC)
# aucs=t(aucs)
scRNASTJT_AUC <- CreateAssayObject(counts = aucs)
scRNASTJT@assays$AUC <- scRNASTJT_AUC
scRNASTJT@assays$AUC@key <- "rna_"
DefaultAssay(scRNASTJT) <- "AUC"
scRNASTJT@meta.data$orig.ident=factor(scRNASTJT@meta.data$orig.ident,levels =c('ST','JT'))
# 画scRNA通路
Idents(scRNASTJT)='orig.ident'
pdf(file="scRNASTJT-Activepathway.pdf",height = 6, width = 12)
DotPlot(scRNASTJT,group.by = 'orig.ident', features = unique(features)) + coord_flip()+labs(x="Pathways")
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()


library(limma)
data=scRNASTJT@assays$AUC@data
dd_CM<- as.data.frame(as.matrix(data))
# colnames(dd_CM)
table(Idents(scRNASTJT))
group2<-  factor(c(rep('JT', 10743),rep('ST', 9321)))#顺序靠后的 vs 顺序靠前的
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
write.table(lrt, 'STvsJT_1.txt', sep = '\t', col.names = NA, quote = FALSE)

GOBP-CELLULAR-RESPONSE-TO-LEPTIN-STIMULUS
GOBP-RESPONSE-TO-LEPTIN
GOBP-NEGATIVE-REGULATION-OF-BROWN-FAT-CELL-DIFFERENTIATION

features=c("WANG-ADIPOGENIC-GENES-REPRESSED-BY-SIRT1",'GOBP-POSITIVE-REGULATION-OF-FAT-CELL-PROLIFERATION','GOBP-CELLULAR-RESPONSE-TO-LEPTIN-STIMULUS','GOBP-RESPONSE-TO-LEPTIN','GOBP-NEGATIVE-REGULATION-OF-BROWN-FAT-CELL-DIFFERENTIATION')
Idents(scRNASTJT)='orig.ident'
scRNASTJT@meta.data$orig.ident=factor(scRNASTJT@meta.data$orig.ident,levels =c('JT','ST'))
pdf(file="scRNASTJT-ActivepathwayFCTOP5.pdf",height = 4, width = 10)
DotPlot(scRNASTJT,group.by = 'orig.ident', features = unique(features)) + coord_flip()+labs(x="Pathways")
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()

WANG-CLASSIC-ADIPOGENIC-TARGETS-OF-PPARG
VERNOCHET-ADIPOGENESIS
BIOCARTA-VOBESITY-PATHWAY
GOBP-NEGATIVE-REGULATION-OF-BROWN-FAT-CELL-DIFFERENTIATION
GOBP-POSITIVE-REGULATION-OF-FAT-CELL-PROLIFERATION

features=c("WANG-CLASSIC-ADIPOGENIC-TARGETS-OF-PPARG",'GOBP-POSITIVE-REGULATION-OF-FAT-CELL-PROLIFERATION','VERNOCHET-ADIPOGENESIS','BIOCARTA-VOBESITY-PATHWAY','GOBP-NEGATIVE-REGULATION-OF-BROWN-FAT-CELL-DIFFERENTIATION')
Idents(scRNASiJi)='orig.ident'
scRNASiJi@meta.data$orig.ident=factor(scRNASiJi@meta.data$orig.ident,levels =c('Ji','Si'))
pdf(file="scRNASiJi-ActivepathwayFCTOP5.pdf",height = 3, width = 10)
DotPlot(scRNASiJi,group.by = 'orig.ident', features = unique(features)) + coord_flip()+labs(x="Pathways")
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()

BIOCARTA-VOBESITY-PATHWAY
WANG-CLASSIC-ADIPOGENIC-TARGETS-OF-PPARG
VERNOCHET-ADIPOGENESIS
GOBP-FAT-PAD-DEVELOPMENT
GOBP-NEGATIVE-REGULATION-OF-BROWN-FAT-CELL-DIFFERENTIATION

features=c("WANG-CLASSIC-ADIPOGENIC-TARGETS-OF-PPARG",'GOBP-FAT-PAD-DEVELOPMENT','VERNOCHET-ADIPOGENESIS','BIOCARTA-VOBESITY-PATHWAY','GOBP-NEGATIVE-REGULATION-OF-BROWN-FAT-CELL-DIFFERENTIATION')
Idents(scRNASeJe)='orig.ident'
scRNASeJe@meta.data$orig.ident=factor(scRNASeJe@meta.data$orig.ident,levels =c('Je','Se'))
pdf(file="scRNASeJe-ActivepathwayFCTOP5.pdf",height = 3, width = 10)
DotPlot(scRNASeJe,group.by = 'orig.ident', features = unique(features)) + coord_flip()+labs(x="Pathways")
# print(DotPlot(obj,group.by = 'cellType', features = unique(features),cluster.idents = T) + coord_flip()+RotatedAxis())
dev.off()
