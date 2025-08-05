##############安装
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("progeny")
## To install the new version until it is submitted to Bioconductor use:
devtools::install_github("saezlab/progeny")
###############加载
library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
library(patchwork)
setwd("C:/Users/LXF/Desktop/fh/copykat")
data<- read.table("strbrafdata674.txt")
gatk<- read.table("gatk.txt")
gatk[1:172,2]<- 'd0'
gatk[173:327,2]<- 'd4'
gatk[328:526,2]<- 'd28'
gatk[527:674,2]<- 'd57'
colnames(gatk)<- c('celltype','stage')
SeuratObject <- CreateSeuratObject(counts = data, project = "mel", meta.data =gatk)
SeuratObject[["percent.mt"]] <- PercentageFeatureSet(SeuratObject, pattern = "^MT-")
SeuratObject <- NormalizeData(SeuratObject)
SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(SeuratObject)
SeuratObject <- ScaleData(SeuratObject, features = all.genes)
SeuratObject <- RunPCA(SeuratObject, features = VariableFeatures(object = SeuratObject))
SeuratObject <- FindNeighbors(SeuratObject,dims = 1:19)
SeuratObject <- FindClusters(SeuratObject, resolution =0.2)#4个，1.1为两个
head(Idents(SeuratObject), 5)
SeuratObject <- RunUMAP(SeuratObject, dims = 1:19)
SeuratObject

###################细胞群通路活性
## Finding differentially expressed features (cluster biomarkers)
SOBJ.markers <- FindAllMarkers(SeuratObject, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25, verbose = FALSE)


## We create a data frame with the specification of the cells that belong to 
## each cluster to match with the Progeny scores.
CellsClusters <- data.frame(Cell = names(Idents(SeuratObject)), 
                            CellType = as.character(Idents(SeuratObject)),
                            stringsAsFactors = FALSE)
DimPlot(SeuratObject, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#####################PROGENy pathway activity scores
## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
# SeuratObject <- progeny(SeuratObject, scale=FALSE, organism="Human", top=500, perm=1, 
#                 return_assay = TRUE)
SeuratObject <- progeny(SeuratObject, scale=FALSE, organism="Mouse", top=500, perm=1, 
                return_assay = TRUE)                

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
SeuratObject <- Seurat::ScaleData(SeuratObject, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(SeuratObject, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
#########plot the different pathway activities for the different cell populations
## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)

# ------------------------------------------------------------------------------------------------


CellsClusters <- data.frame(Cell = names(Idents(scRNA)), 
                            CellType = as.character(Idents(scRNA)),
                            stringsAsFactors = FALSE)
DimPlot(scRNA, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#####################PROGENy pathway activity scores
## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
# scRNA <- progeny(scRNA, scale=FALSE, organism="Human", top=500, perm=1, 
#                 return_assay = TRUE)
scRNA <- progeny(scRNA, scale=FALSE, organism="Mouse", top=500, perm=1, 
                return_assay = TRUE)                

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
scRNA <- Seurat::ScaleData(scRNA, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(scRNA, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
#########plot the different pathway activities for the different cell populations
## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)

ggsave('progeny_hmap.pdf',progeny_hmap,width=10,height=10)



DefaultAssay(scRNA)<-"progeny"
scRNASeJeprogeny=subset(x=scRNA,subset=(orig.ident == "Se"|orig.ident == "Je"))
p1= FeaturePlot(scRNASeJeprogeny,features ="VEGF",coord.fixed=T,order=T,cols=c('#6a5acd','#cb5cbd'),split.by='orig.ident')
ggsave('scRNASeJeprogeny.pdf',p1)



# -------------------------------


scRNASeJe=subset(x=scRNA,subset=(orig.ident == "Se"|orig.ident == "Je"))
scRNASeJeEC=subset(x=scRNASeJe,subset=(cellType == "EC"))
Idents(scRNASeJeEC)='orig.ident'
CellsClusters <- data.frame(Cell = names(Idents(scRNASeJeEC)), 
                            CellType = as.character(Idents(scRNASeJeEC)),
                            stringsAsFactors = FALSE)
# DimPlot(scRNASeJeEC, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#####################PROGENy pathway activity scores
## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
# scRNASeJeEC <- progeny(scRNASeJeEC, scale=FALSE, organism="Human", top=500, perm=1, 
#                 return_assay = TRUE)
scRNASeJeEC <- progeny(scRNASeJeEC, scale=FALSE, organism="Mouse", top=500, perm=1, 
                return_assay = TRUE)                

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
scRNASeJeEC <- Seurat::ScaleData(scRNASeJeEC, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(scRNASeJeEC, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
#########plot the different pathway activities for the different cell populations
## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df),fontsize=20, 
                        fontsize_row = 20, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy", angle_col = 0,
                        treeheight_col = 0,  border_color = NA)

ggsave('SeJe-EC-progeny_hmap.pdf',progeny_hmap,width=4,height=6)

# GSVA
genelist<- list.files('./gmtps')

#Mm_m5 <- msigdbr(species = "Mus musculus", category = "m5")## 定义基因集，选取C2
#fgsea_sets = mdb_c2 [grep("^GO",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets<- list()
for (i in 1:length(genelist)) {
  dir<- paste0('./gmtps/',genelist[i])
  temp<- read.gmt(dir)
  fgsea_sets[[names(table(temp$term))]]<- temp$gene
  
}
scRNASeJemyoFB=subset(x=scRNASeJe,subset=(cellType == "myoFB"))
data=scRNASeJemyoFB@assays$RNA@data

expr=as.matrix(data) 
kegg <- gsva(expr, fgsea_sets, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
write.csv(file='kegg.csv',kegg)
# p=pheatmap(kegg)#绘制热图
Idents(scRNASeJemyoFB)='orig.ident'
CellsClusters <- data.frame(Cell = names(Idents(scRNASeJemyoFB)), 
                            CellType = as.character(Idents(scRNASeJemyoFB)),
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
ggsave('SeJe-myoFB-gsva_hmap.pdf',progeny_hmap,width=10,height=6)

progeny_hmap = pheatmap(t(summarized_kegg_scores_df),fontsize=16, 
                        fontsize_row = 16, 
                        color=myColor, 
                        main = "GSVA", angle_col = 0,
                        treeheight_col = 0,  border_color = NA)
ggsave('SeJe-myoFB-gsva_hmap2.pdf',progeny_hmap,width=10,height=6)

# GSVA
genelist<- list.files('./gmtps')

#Mm_m5 <- msigdbr(species = "Mus musculus", category = "m5")## 定义基因集，选取C2
#fgsea_sets = mdb_c2 [grep("^GO",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets<- list()
for (i in 1:length(genelist)) {
  dir<- paste0('./gmtps/',genelist[i])
  temp<- read.gmt(dir)
  fgsea_sets[[names(table(temp$term))]]<- temp$gene
  
}
scRNASeJeEC=subset(x=scRNASeJe,subset=(cellType == "EC"))
data=scRNASeJeEC@assays$RNA@data

expr=as.matrix(data) 
kegg <- gsva(expr, fgsea_sets, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
write.csv(file='kegg.csv',kegg)
# p=pheatmap(kegg)#绘制热图
Idents(scRNASeJeEC)='orig.ident'
CellsClusters <- data.frame(Cell = names(Idents(scRNASeJeEC)), 
                            CellType = as.character(Idents(scRNASeJeEC)),
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
ggsave('SeJe-EC-gsva_hmap.pdf',progeny_hmap,width=10,height=6)

progeny_hmap = pheatmap(t(summarized_kegg_scores_df),fontsize=16, 
                        fontsize_row = 16, 
                        color=myColor,  
                        main = "GSVA", angle_col = 0,
                        treeheight_col = 0,  border_color = NA)
ggsave('SeJe-EC-gsva_hmap2.pdf',progeny_hmap,width=10,height=6)


df=data.frame(name=c('EC_VEGF','EC_VEGF','EC_VEGF','EC_VEGF'),Type=c('Nucb2_s','Vegfa_s','DoubleT','DoubleF'),Num=c(111,23,15,598))
data3 = ddply(df,'name',transform,percent_Num=Num/sum(Num)*100)
p <- ggplot(data3,aes(x=name,y=percent_Num,fill=Type))+geom_bar(stat = 'identity',width = 0.5,colour='black')+theme_classic()+labs(x='',y='Percentage')+theme(axis.title = element_text(size=20),axis.text = element_text(size=20))+scale_y_continuous(breaks=seq(0,100,25),labels=c('0','25%','50%','75%','100%'))+geom_text(aes(label=paste0(sprintf("%.1f", percent_Num), "%")),position = position_stack(vjust = 0.5), size = 8,color="black")+scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72"))+theme( legend.text = element_text(size = 18))
p

df=data.frame(name=c('myoFB_IGF','myoFB_IGF','myoFB_IGF','myoFB_IGF'),Type=c('Nucb2_s','igf1r_s','DoubleT','DoubleF'),Num=c(378,1874,2144,421))
data3 = ddply(df,'name',transform,percent_Num=Num/sum(Num)*100)
p <- ggplot(data3,aes(x=name,y=percent_Num,fill=Type))+geom_bar(stat = 'identity',width = 0.5,colour='black')+theme_classic()+labs(x='',y='Percentage')+theme(axis.title = element_text(size=20),axis.text = element_text(size=20))+scale_y_continuous(breaks=seq(0,100,25),labels=c('0','25%','50%','75%','100%'))+geom_text(aes(label=paste0(sprintf("%.1f", percent_Num), "%")),position = position_stack(vjust = 0.5), size = 8,color="black")+scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72"))+theme( legend.text = element_text(size = 18))
p

df=data.frame(name=factor(c(rep('Se',8),rep('Je',8))),Type=c('0','1','2','3','4','5','6','7','0','1','2','3','4','5','6','7'),Num=c(6,27,403,117,8,2,49,0,1578,1440,407,588,57,50,2,40))
data3 = ddply(df,'name',transform,percent_Num=Num/sum(Num)*100)
p <- ggplot(data3,aes(x=name,y=percent_Num,fill=Type))+geom_bar(stat = 'identity',width = 0.5,colour='black')+theme_classic()+labs(x='',y='Percentage')+theme(axis.title = element_text(size=20),axis.text = element_text(size=20))+scale_y_continuous(breaks=seq(0,100,25),labels=c('0','25%','50%','75%','100%'))+geom_text(aes(label=paste0(sprintf("%.1f", percent_Num), "%")),position = position_stack(vjust = 0.5), size = 8,color="black")+scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#9932CC","#8B008B","#8B4513","#DEB887"))+theme( legend.text = element_text(size = 18))
p

df=data.frame(name=factor(c(rep('Se',8),rep('Je',8))),Type=c('0','1','2','3','4','5','6','7','0','1','2','3','4','5','6','7'),Num=c(6,27,403,117,8,2,49,0,1578,1440,407,588,57,50,2,40))
data3 = ddply(df,'name',transform,percent_Num=Num/sum(Num)*100)
data3$name=factor(data3$name,levels=c('Se','Je'))
write.csv(data3,file='sub_data.csv')
p <- ggplot(data3,aes(x=name,y=percent_Num,fill=Type))+geom_bar(stat = 'identity',width = 0.5,colour='black')+theme_classic()+labs(x='',y='Percentage')+theme(axis.title = element_text(size=20),axis.text = element_text(size=20))+scale_y_continuous(breaks=seq(0,100,25),labels=c('0','25%','50%','75%','100%'))+scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#9932CC","#8B008B","#8B4513","#DEB887"))+theme( legend.text = element_text(size = 18))
p
