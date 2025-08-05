#devtools::install_github('sqjin/CellChat')
# devtools::install_github('jinworks/CellChat')
# devtools::install_github('immunogenomics/presto')
library(presto)
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)
TNK$new<- paste0(TNK$response,'_',TNK$stage)
table(TNK$new)
norespost<- subset(TNK,subset=new=='no-res_post')
norespre<- subset(TNK,subset=new=='no-res_pre')
respost<- subset(TNK,subset=new=='res_post')
respre <- subset(TNK,subset=new=='res_pre')

"CLZBMMC531" "CLZBMMC620" "CLZPBMC531" "CLZPBMC620" "LBHBMMC616" 
"LBHBMMC702" "LBHPBMC616" "LBHPBMC702" "MHJBMMC711"

data.input = respre@assays$RNA@data # normalized data matrix
meta = respre@meta.data # a dataframe with rownames containing cell mata data


###创建CellChat 对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellType")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cellType") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 10) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

###细胞通信网络的推断
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
cellchat <- updateCellChat(cellchat)
table(cellchat@idents)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Interaction weights/strength")


save(cellchat,file = 'D:\\wmy\\project\\ALL\\cellchat2\\cell-cell\\LBHBMMC616.RData')

"""
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]  
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
"""
##WEIGHT
mat <- cellchat@net$weight
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
rownames(mat)
i=31
mat2[i, ] <- mat[i, ] 
mat2[,i] <- mat[,i] 
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])

##COUNT
mat <- cellchat@net$count
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
rownames(mat)
i=31
mat2[i, ] <- mat[i, ] 
mat2[,i] <- mat[,i] 
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])

cellchat@meta$new3<- factor(cellchat@meta$new3,
                            levels = c("B","Fibroblast","progenitor","Macrophage","pDC","DN","DPbla-1",
                                       "DPbla-2","DPbla-3","DPbla-4","DPbla-5","DPbla-6","DPbla-7",
                                       "DPre-1","DPre-2","DPre-3","DPre-4","DPre-5","DPre-6","DPre-7",
                                       "DPre-8","DPre-9","DPre-10","DPSel-1","DPSel-2","DPSel-3","DPSel-4",
                                       "DPSel-5","DPSel-6","DPSel-7","CD4 SP-1","CD4 SP-2","CD4 SP-3",
                                       "CD4 SP-4","CD4 Treg","CD8 SP-1","CD8 SP-2"))



###细胞通信网络的可视化
cellchat@netP$pathways
[1] "MHC-I"       "MHC-II"      "PECAM1"      "LCK"         "SELPLG"      "ITGAL-ITGB2" "CD86"        "CD22"       
[9] "CD45"        "ICAM"        "CD52"        "ALCAM"       "CD6"         "CD226"       "PVR"         "THY1"       
[17] "CD96"        "SEMA4"       "APP"         "CD80"        "SELL"        "ICOS"        "VCAM"        "PD-L1"      
[25] "CADM"        "CD48"        "NECTIN"      "NOTCH"       "CD23"        "CLEC"        "NKG2D"       "CDH1"       
[33] "CD200"       "EPHA"        "BST2"        "CDH"         "JAM"         "PTPRM"       "PDL2"        "CDH5"       
[41] "EPHB"        "LAIR1"       "CD40"        "SEMA5"       "SEMA6"       "CD39"        "CD46"        "NEGR"       
[49] "MPZ"
pathways.show <- c("MHC-I") # Hierarchy plot# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 

levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents))
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)


netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")#> Do heatmap based on a single ohttp://127.0.0.1:18459/graphics/plot_zoom_png?width=1200&height=900bject

# Chord diagram
group.cellType <- c("DPSel-7","pDC") 
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))#> Plot the aggregated cell-cell communication network at the signaling pathway level#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

#> [[1]]# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# Access all the signaling pathways showing significant communicationspathways.show.all <- cellchat@netP$pathways# check the order of cell identity to set suitable 
vertex.receiverlevels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {  
  # Visualize communication network associated with both signaling pathway and individual L-R pairs  
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway  
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])  
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
  }

###可视化由多个配体受体或信号通路调节的细胞通信
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
fic1<-netVisual_bubble(cellchat, sources.use = 31, targets.use = c(1:38), remove.isolate = FALSE)#> Comparing communications on a single object
fic2<-netVisual_bubble(cellchat, sources.use =c(1:38), targets.use = 31, remove.isolate = FALSE)#
fic1/fic2
###和弦图
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 31, targets.use = c(1:4), lab.cex = 0.5,legend.pos.y = 30)#> Note: The first link end is drawn out of sector 'MIF'.
# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = 31, targets.use = c(1:4), legend.pos.x = 15)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = 31, targets.use = c(1:4), signaling = c("CCL","CXCL"),legend.pos.x = 8)#> Note: The second link end is drawn out of sector 'CXCR4 '.#> Note: The first link end is drawn out of sector 'CXCL12 '.
# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)#> Note: The second link end is drawn out of sector ' '.#> Note: The first link end is drawn out of sector 'MIF'.#> Note: The second link end is drawn out of sector ' '.#> Note: The first link end is drawn out of sector 'CXCL '.

###使用小提琴/点图绘制信号基因表达分布
plotGeneExpression(cellchat, signaling = "MHC-I")#> Registered S3 method overwritten by 'spatstat':#>   method     from#>   print.boxx cli#> Scale for 'y' is already present. Adding another scale for 'y', which will#> replace the existing scale.#> Scale for 'y' is already present. Adding another scale for 'y', which will#> replace the existing scale.#> Scale for 'y' is already present. Adding another scale for 'y', which will#> replace the existing scale.
plotGeneExpression(cellchat, signaling = "MHC-I", enriched.only = FALSE)


###细胞通信网络系统分析
# Compute the network centrality 
scorescellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(scorescellchat, signaling = pathways.show, width = 16, height = 12, font.size = 10)

###在 2D 空间中可视化占主导地位的发送器（源）和接收器（目标）
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(scorescellchat)#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(scorescellchat, signaling = c("CXCL", "CCL"))#> Signaling role analysis on the cell-cell communication network from user's inputgg1 + gg2
gg1+gg2

###识别对某些细胞组的传出或传入信号贡献最大的信号
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(scorescellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(scorescellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(scorescellchat, signaling = c( "MHC-I","MHC-II" ))


###确定全局通信模式，探索多个细胞类型和信号通路如何协调在一起
library(NMF)#> Loading required package: pkgmaker#> Loading required package: registry#> Loading required package: rngtools#> Loading required package: cluster#> NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 15/16#>   To enable shared memory capabilities, try: install.extras('#> NMF#> ')#> #> Attaching package: 'NMF'#> The following objects are masked from 'package:igraph':#> #>     algorithm, comparelibrary(ggalluvial)
selectK(cellchat, pattern = "outgoing")

nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
library(ggalluvial)
netAnalysis_river(cellchat, pattern = "outgoing")#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

###识别和可视化目标细胞的传入通信模式
selectK(cellchat, pattern = "incoming")

nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

###信号网络的多重和分类学习分析
cellchat <- computeNetSimilarity(cellchat, type = "functional")
reticulate::py_install(packages = 'umap-learn')
cellchat <- netEmbedding(cellchat, type = "functional")#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")#> Classification learning of the signaling networks for a single dataset# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

################以上为历史代码

####细胞通讯####
#devtools::install_github('sqjin/CellChat')
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)
# 读取数据
scRNA=readRDS('scRNA_seurat_cellType.RDS')
# 拆分样本
scRNASe=subset(x = scRNA, subset = (orig.ident == "Se"))
scRNAJe=subset(x = scRNA, subset = (orig.ident == "Je"))
scRNASi=subset(x = scRNA, subset = (orig.ident == "Si"))
scRNAJi=subset(x = scRNA, subset = (orig.ident == "Ji"))
scRNAST=subset(x = scRNA, subset = (orig.ident == "ST"))
scRNAJT=subset(x = scRNA, subset = (orig.ident == "JT"))

# CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data 
# showDatabaseCategory(CellChatDB)
i=1
samplenames=c('Se','Si','ST','Je','Ji','JT')
searchs=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
for (seuratobj_raw in c(scRNASe,scRNASi,scRNAST,scRNAJe,scRNAJi,scRNAJT)){
  seuratobj <- seuratobj_raw
  samplename=samplenames[i]
  print(samplename)
  print(seuratobj)
  # i=i+1
  data.input = seuratobj@assays$RNA@data # normalized data matrix
  meta = seuratobj@meta.data # a dataframe with rownames containing cell mata data
  ###创建CellChat 对象
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellType")
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "cellType") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  # CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data 
  CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data 
  showDatabaseCategory(CellChatDB)

  # Show the structure of the database
  dplyr::glimpse(CellChatDB$interaction)

  for (search in searchs){
    # use a subset of CellChatDB for cell-cell communication analysis
    CellChatDB.use <- subsetDB(CellChatDB, search = search) # use Secreted Signaling
    # use all CellChatDB for cell-cell communication analysis
    # CellChatDB.use <- CellChatDB # simply use the default CellChatDB

    # set the used database in the object
    cellchat@DB <- CellChatDB.use

    cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
    future::plan("multisession", workers = 50) # do parallel

    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.mouse)

    ###细胞通信网络的推断
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)

    groupSize <- as.numeric(table(cellchat@idents))

    pdf(file=sprintf("%s%s_netVisual_circle.pdf",samplename,search),height = 15, width = 15)
    par(mfrow = c(1,2), xpd=TRUE)
    cellchat <- updateCellChat(cellchat)
    table(cellchat@idents)
    
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                    weight.scale = T, label.edge= F, title.name = "Number of interactions")

    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Interaction strength")
    dev.off()
 
    save(cellchat,file = sprintf('%s_cellchat_%s.RData',samplename,search))}

  i=i+1
  }

for (samplename in samplenames){
  for (search in searchs){
    load(sprintf('%s_cellchat_%s.RData',samplename,search))
    # 重新绘图
    pdf(file=sprintf("%s%s_netVisual_circle.pdf",samplename,search),height = 15, width = 15)
    par(mfrow = c(1,2), cex.main = 2, xpd=TRUE)
    
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                    weight.scale = T, label.edge= F, title.name = "Number of interactions")

    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Interaction weights/strength")
    dev.off()

  }
}

# # 整合所有cellchat
# object.list <- list(Se_cellchat=Se_cellchat,Si_cellchat=Si_cellchat, ST_cellchat = ST_cellchat,Je_cellchat = Je_cellchat,Ji_cellchatt=Ji_cellchat,JT_cellchat=JT_cellchat)
# #run netAnalysis_computeCentrality
# object.list<- lapply(object.list,function(x){
#   x=netAnalysis_computeCentrality(x)
# })
# cellchat_m <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)

# ###取子集
# names(object.list)
# cellchat_m <- setIdent(cellchat_m, ident.use = "sub_clusters")#以sub_clusters设置细胞标志
# table(cellchat_m@idents)
# cellchat_m
# cellchat_p<- subsetCellChat(cellchat_m,idents.use =c('DPSel-7','pDC'))#根据sub_clusters中的细胞标志提取'DPSel-7','pDC'这两种类型数据

# 整合Se Je Secreted Signaling cellchat
load('Se_cellchat_Secreted Signaling.RData')
Se_cellchat=cellchat
load('Je_cellchat_Secreted Signaling.RData')
Je_cellchat=cellchat
object.list <- list(Se_cellchat=Se_cellchat,Je_cellchat = Je_cellchat)
#run netAnalysis_computeCentrality
object.list<- lapply(object.list,function(x){
  x=netAnalysis_computeCentrality(x)
})
cellchat_m <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)

######part1 从宏观角度预测细胞通讯
###比较交互总数和交互强度
pdf(file="SeJeChatNumall.pdf",height = 8, width = 8)
compareInteractions(cellchat_m, show.legend = F, group = c(1,2),size.text = 20)#group为相互比较的不同的cellchat，通过names(object.list)查看用来比较的cellchat所在的位置
compareInteractions(cellchat_m, show.legend = F, group = c(1,2), measure = "weight",size.text = 20)
dev.off()
###不同细胞群之间的相互作用数量或强度的差异 circle
pdf(file="SeJeChatNumCell.pdf",height = 8, width = 8)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = c(1,2))
netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = c(1,2))
dev.off()
###不同细胞群之间的相互作用数量或强度的差异 heatmap
pdf(file="SeJeChatNumCellH.pdf",height = 8, width = 8)
netVisual_heatmap(cellchat_m,comparison = c(1,2), font.size = 20, font.size.title = 20)
#> Do heatmap based on a merged object
netVisual_heatmap(cellchat_m, measure = "weight",comparison = c(1,2), font.size = 20, font.size.title = 20)
dev.off()


library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

scRNA=readRDS('scRNA_seurat_cellType.RDS')
scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))
scRNASeJemyoFB=subset(x=scRNASeJe,subset=(cellType == "myoFB"))
scRNASeJeEC=subset(x=scRNASeJe,subset=(cellType == "EC"))
scRNASeJeAP=subset(x=scRNASeJe,subset=(cellType == "AP"))
# AP亚群分类
DefaultAssay(scRNASeJeAP) <- "RNA"

scRNASeJeAP <- NormalizeData(scRNASeJeAP, normalization.method = "LogNormalize", scale.factor = 10000)

scRNASeJeAP <- FindVariableFeatures(scRNASeJeAP,nfeatures = 3000)
scRNASeJeAP <- ScaleData(scRNASeJeAP,vars.to.regress = 'percent.mt')
scRNASeJeAP <- RunPCA(scRNASeJeAP)
scRNASeJeAP <- FindNeighbors(scRNASeJeAP, reduction = "pca", dims = 1:30)
scRNASeJeAP <- FindClusters(scRNASeJeAP, resolution = 0.3)
scRNASeJeAP <- RunUMAP(scRNASeJeAP, reduction = "pca", dims = 1:30) 
library(harmony)
scRNASeJeAP <- RunHarmony(scRNASeJeAP, group.by.vars = "orig.ident")
scRNASeJeAP  <- FindNeighbors(scRNASeJeAP , dims = 1:30 , reduction = "harmony")
scRNASeJeAP  <- FindClusters(scRNASeJeAP, save.snn=T , resolution = 0.3)
scRNASeJeAP  <- RunUMAP(scRNASeJeAP , dims=1:30,reduction='harmony')
saveRDS(scRNASeJeAP,'scRNASeJeAPharmony.rda')
# scRNASeJeAP=readRDS('scRNASeJeAPharmony.rda')
###AP亚群细胞通讯
# scRNASeJeAP6=subset(x=scRNASeJeAP,subset=(seurat_clusters== "6"))
scRNASeJeAP$new='Nucb2_0'
# filtered_cells <- WhichCells(scRNASeJeAP, expression = Nucb2 == 0)
Idents(scRNASeJeAP)='new'
filtered_cells <- WhichCells(scRNASeJeAP, expression = Nucb2 > 0)
Idents(scRNASeJeAP, cells = filtered_cells) <- 'Nucb2_1'
scRNASeJeAP@meta.data$new=factor(Idents(scRNASeJeAP))
table(Idents(scRNASeJeAP))


# seuratobj=scRNASeJeAP
samplename='scRNASeJeAP'
search='Secreted Signaling'

data.input = scRNASeJeAP@assays$RNA@data # normalized data matrix
meta = scRNASeJeAP@meta.data # a dataframe with rownames containing cell mata data
###创建CellChat 对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "new")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "new") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
# CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data 
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data 
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = search) # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 50) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj=PPI.mouse)

###细胞通信网络的推断
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

save(cellchat,file = sprintf('%s_cellchat_%s.RData',samplename,search))




scRNASeJeAPmyoFBEC=subset(x=scRNASeJe,subset=(cellType == "AP"|cellType == "myoFB"|cellType == "EC"))

Idents(scRNASeJeAP)='new'
Idents(scRNASeJeAPmyoFBEC)='cellType'

filtered_cells <- WhichCells(scRNASeJeAP, expression = Nucb2 > 0)
# filtered_cells <- WhichCells(scRNASeJeAP, expression = Nucb2 > 0&seurat_clusters='6')
Idents(scRNASeJeAPmyoFBEC, cells = filtered_cells) <- 'APNucb2_1'


filtered_cells <- WhichCells(scRNASeJeAP, expression = Nucb2 == 0)
Idents(scRNASeJeAPmyoFBEC, cells = filtered_cells) <- 'APNucb2_0'

scRNASeJeAPmyoFBEC@meta.data$cellType=factor(Idents(scRNASeJeAPmyoFBEC))
table(Idents(scRNASeJeAPmyoFBEC))

scRNASeAPmyoFBEC=subset(x = scRNASeJeAPmyoFBEC, subset = (orig.ident == "Se"))

samplename='scRNASeAPmyoFBEC'
search='Secreted Signaling'

data.input = scRNASeAPmyoFBEC@assays$RNA@data # normalized data matrix
meta = scRNASeAPmyoFBEC@meta.data # a dataframe with rownames containing cell mata data
###创建CellChat 对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellType")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cellType") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
# CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data 
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data 
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = search) # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 50) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj=PPI.mouse)

###细胞通信网络的推断
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

save(cellchat,file = sprintf('%s_cellchat_%s.RData',samplename,search))

scRNAJeAPmyoFBEC=subset(x = scRNASeJeAPmyoFBEC, subset = (orig.ident == "Je"))
samplename='scRNAJeAPmyoFBEC'
search='Secreted Signaling'

data.input = scRNAJeAPmyoFBEC@assays$RNA@data # normalized data matrix
meta = scRNAJeAPmyoFBEC@meta.data # a dataframe with rownames containing cell mata data
###创建CellChat 对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellType")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cellType") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
# CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data 
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data 
# showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = search) # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 50) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj=PPI.mouse)

###细胞通信网络的推断
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

save(cellchat,file = sprintf('%s_cellchat_%s.RData',samplename,search))



#对细胞交互进行展示

load(sprintf('%s_cellchat_%s.RData','scRNASeAPmyoFBEC',search))
Se_cellchat=cellchat
groupSize <- as.numeric(table(cellchat@idents))
pdf(file=sprintf("Se_netVisual_circle_%s.pdf",abbr),height = 8, width = 8)
par(mfrow = c(1,2), xpd=TRUE)
cellchat <- updateCellChat(cellchat)
table(cellchat@idents)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                weight.scale = T, label.edge= F, title.name = "Number of inferred interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Interaction strength")
dev.off()

load(sprintf('%s_cellchat_%s.RData','scRNAJeAPmyoFBEC',search))
Je_cellchat=cellchat
groupSize <- as.numeric(table(cellchat@idents))
pdf(file=sprintf("Je_netVisual_circle_%s.pdf",abbr),height = 8, width = 8)
par(mfrow = c(1,2), xpd=TRUE)
cellchat <- updateCellChat(cellchat)
table(cellchat@idents)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                weight.scale = T, label.edge= F, title.name = "Number of inferred interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Interaction strength")
dev.off()

object.list <- list(Je=Je_cellchat,Se=Se_cellchat)
#run netAnalysis_computeCentrality
object.list<- lapply(object.list,function(x){
    x=netAnalysis_computeCentrality(x)})
cellchat_m <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
####从宏观角度预测细胞通讯####
#比较交互总数和交互强度

compare='subSeJe'
abbr='SS'
group=c(1,2)
pdf(file=sprintf("%sChatNumall_%s.pdf",compare,abbr),height = 4.5, width = 4)
par(mfrow = c(1,1), xpd=TRUE)
compareInteractions(cellchat_m, show.legend = F, group = (1:2),size.text = 20)#group颜色向量 默认measure='count'
# print(compareInteractions(cellchat_m, show.legend = F, group = (1:2),size.text = 20))
compareInteractions(cellchat_m, show.legend = F, group = (1:2), measure = "weight",size.text = 20)
# print(compareInteractions(cellchat_m, show.legend = F, group = (1:2), measure = "weight",size.text = 20))
dev.off()

######不同细胞群之间的相互作用数量或强度的差异 circle
pdf(file=sprintf("%sChatNumCell_%s.pdf",compare,abbr),height = 8, width = 8)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group)#group为相互比较的不同的cellchat，通过names(object.list)查看用来比较的cellchat所在的位置
# print(netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group))
netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group)
# print(netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group))
dev.off()

####识别上调和下调的信号配体对####
pdf(file=sprintf("%sChatbubbleall_%s.pdf",compare,abbr),height = 8, width = 8)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_bubble(cellchat_m,comparison = group, angle.x = 90, font.size = 16) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 16))
dev.off()


########AP6 AP3 myoFB相互作用
scRNA=readRDS('../scRNA_seurat_cellType.RDS')
scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))
scRNASeJeAP=readRDS('scRNASeJeAPharmony.rda')
###AP亚群细胞通讯
scRNASeJeAPmyoFB=subset(x=scRNASeJe,subset=(cellType == "AP"|cellType == "myoFB"|cellType == "EC"))
Idents(scRNASeJeAPmyoFB)='cellType'
filtered_cellsAP6 <- WhichCells(scRNASeJeAP,idents='6')
Idents(scRNASeJeAPmyoFB, cells = filtered_cellsAP6) <- 'AP6'
filtered_cellsAP3 <- WhichCells(scRNASeJeAP, idents ='3')
Idents(scRNASeJeAPmyoFB, cells = filtered_cellsAP3) <- 'AP3'
scRNASeJeAPmyoFB@meta.data$cellType=factor(Idents(scRNASeJeAPmyoFB),levels=c('AP3','AP6','myoFB','EC'))
scRNASeJeAP63myoFB=subset(x=scRNASeJeAPmyoFB,subset=(cellType == "AP6"|cellType == "AP3"|cellType == "myoFB"|cellType == "EC"))

scRNASeAP63myoFB=subset(x = scRNASeJeAP63myoFB, subset = (orig.ident == "Se"))
samplename='scRNASeAP63myoFB'
search='Secreted Signaling'
data.input = scRNASeAP63myoFB@assays$RNA@data # normalized data matrix
meta = scRNASeAP63myoFB@meta.data # a dataframe with rownames containing cell mata data
###创建CellChat 对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellType")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cellType") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data 
# showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = search) # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 50) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj=PPI.mouse)
###细胞通信网络的推断
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 0)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
save(cellchat,file = sprintf('%s_cellchat_%s.RData',samplename,search))

scRNAJeAP63myoFB=subset(x = scRNASeJeAP63myoFB, subset = (orig.ident == "Je"))
samplename='scRNAJeAP63myoFB'
search='Secreted Signaling'
data.input = scRNAJeAP63myoFB@assays$RNA@data # normalized data matrix
meta = scRNAJeAP63myoFB@meta.data # a dataframe with rownames containing cell mata data
###创建CellChat 对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellType")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cellType") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data 
# showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = search) # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 50) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj=PPI.mouse)
###细胞通信网络的推断
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 0)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
save(cellchat,file = sprintf('%s_cellchat_%s.RData',samplename,search))

load(sprintf('%s_cellchat_%s.RData','scRNASeAP63myoFB',search))
Se_cellchat=cellchat
groupSize <- as.numeric(table(cellchat@idents))
pdf(file=sprintf("Se-AP63_netVisual_circle_%s.pdf",abbr),height = 8, width = 8)
par(mfrow = c(1,2), xpd=TRUE)
cellchat <- updateCellChat(cellchat)
table(cellchat@idents)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                weight.scale = T, label.edge= F, title.name = "Number of inferred interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Interaction strength")
dev.off()

load(sprintf('%s_cellchat_%s.RData','scRNAJeAP63myoFB',search))
Je_cellchat=cellchat
groupSize <- as.numeric(table(cellchat@idents))
pdf(file=sprintf("Je-AP63_netVisual_circle_%s.pdf",abbr),height = 8, width = 8)
par(mfrow = c(1,2), xpd=TRUE)
cellchat <- updateCellChat(cellchat)
table(cellchat@idents)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                weight.scale = T, label.edge= F, title.name = "Number of inferred interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Interaction strength")
dev.off()

object.list <- list(Je=Je_cellchat,Se=Se_cellchat)

#run netAnalysis_computeCentrality
object.list<- lapply(object.list,function(x){
    x=netAnalysis_computeCentrality(x)})

cellchat_m <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
####从宏观角度预测细胞通讯####
#比较交互总数和交互强度

compare='subSeJeAP63'
abbr='SS'
group=c(1,2)
pdf(file=sprintf("%sChatNumall_%s.pdf",compare,abbr),height = 4.5, width = 4)
par(mfrow = c(1,1), xpd=TRUE)
compareInteractions(cellchat_m, show.legend = F, group = (1:2),size.text = 20)#group颜色向量 默认measure='count'
# print(compareInteractions(cellchat_m, show.legend = F, group = (1:2),size.text = 20))
compareInteractions(cellchat_m, show.legend = F, group = (1:2), measure = "weight",size.text = 20)
# print(compareInteractions(cellchat_m, show.legend = F, group = (1:2), measure = "weight",size.text = 20))
dev.off()

######不同细胞群之间的相互作用数量或强度的差异 circle
pdf(file=sprintf("%sChatNumCell_%s.pdf",compare,abbr),height = 4, width = 4)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group,sources.use=c('AP3','AP6','myoFB','EC'),targets.use=c('AP3','AP6','myoFB','EC'))#group为相互比较的不同的cellchat，通过names(object.list)查看用来比较的cellchat所在的位置
# print(netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group))
netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group,sources.use=c('AP3','AP6','myoFB','EC'),targets.use=c('AP3','AP6','myoFB','EC'))
# print(netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group))
dev.off()

####识别上调和下调的信号配体对####
pdf(file=sprintf("%sChatbubbleall_%s.pdf",compare,abbr),height = 7, width = 10)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_bubble(cellchat_m,comparison = group,sources.use=c('AP3','AP6','myoFB','EC'),targets.use=c('AP3','AP6','myoFB','EC'), angle.x = 90, font.size = 16) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 16))
dev.off()


########AP6 AP3 EC相互作用
scRNA=readRDS('../scRNA_seurat_cellType.RDS')
scRNASeJe=subset(x = scRNA, subset = (orig.ident == "Se"|orig.ident == "Je"))
scRNASeJeAP=readRDS('scRNASeJeAPharmony.rda')
###AP亚群细胞通讯
scRNASeJeAPEC=subset(x=scRNASeJe,subset=(cellType == "AP"|cellType == "EC"))
Idents(scRNASeJeAPEC)='cellType'
filtered_cellsAP6 <- WhichCells(scRNASeJeAP,idents='6')
Idents(scRNASeJeAPEC, cells = filtered_cellsAP6) <- 'AP6'
filtered_cellsAP3 <- WhichCells(scRNASeJeAP, idents ='3')
Idents(scRNASeJeAPEC, cells = filtered_cellsAP3) <- 'AP3'
scRNASeJeAPEC@meta.data$cellType=factor(Idents(scRNASeJeAPEC),levels=c('AP3','AP6','EC'))
scRNASeJeAP63EC=subset(x=scRNASeJeAPEC,subset=(cellType == "AP6"|cellType == "AP3"|cellType == "EC"))

scRNASeAP63EC=subset(x = scRNASeJeAP63EC, subset = (orig.ident == "Se"))
samplename='scRNASeAP63EC'
search='Secreted Signaling'
data.input = scRNASeAP63EC@assays$RNA@data # normalized data matrix
meta = scRNASeAP63EC@meta.data # a dataframe with rownames containing cell mata data
###创建CellChat 对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellType")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cellType") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data 
# showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = search) # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 50) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj=PPI.mouse)
###细胞通信网络的推断
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 0)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
save(cellchat,file = sprintf('%s_cellchat_%s.RData',samplename,search))

scRNAJeAP63EC=subset(x = scRNASeJeAP63EC, subset = (orig.ident == "Je"))
samplename='scRNAJeAP63EC'
search='Secreted Signaling'
data.input = scRNAJeAP63EC@assays$RNA@data # normalized data matrix
meta = scRNAJeAP63EC@meta.data # a dataframe with rownames containing cell mata data
###创建CellChat 对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellType")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cellType") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data 
# showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = search) # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 50) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj=PPI.mouse)
###细胞通信网络的推断
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 0)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
save(cellchat,file = sprintf('%s_cellchat_%s.RData',samplename,search))

load(sprintf('%s_cellchat_%s.RData','scRNASeAP63EC',search))
Se_cellchat=cellchat
groupSize <- as.numeric(table(cellchat@idents))
pdf(file=sprintf("Se-AP63EC_netVisual_circle_%s.pdf",abbr),height = 8, width = 8)
par(mfrow = c(1,2), xpd=TRUE)
cellchat <- updateCellChat(cellchat)
table(cellchat@idents)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                weight.scale = T, label.edge= F, title.name = "Number of inferred interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Interaction strength")
dev.off()

load(sprintf('%s_cellchat_%s.RData','scRNAJeAP63EC',search))
Je_cellchat=cellchat
groupSize <- as.numeric(table(cellchat@idents))
pdf(file=sprintf("Je-AP63EC_netVisual_circle_%s.pdf",abbr),height = 8, width = 8)
par(mfrow = c(1,2), xpd=TRUE)
cellchat <- updateCellChat(cellchat)
table(cellchat@idents)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                weight.scale = T, label.edge= F, title.name = "Number of inferred interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Interaction strength")
dev.off()

object.list <- list(Je=Je_cellchat,Se=Se_cellchat)

#run netAnalysis_computeCentrality
object.list<- lapply(object.list,function(x){
    x=netAnalysis_computeCentrality(x)})

cellchat_m <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
####从宏观角度预测细胞通讯####
#比较交互总数和交互强度

compare='subSeJeAP63EC'
abbr='SS'
group=c(1,2)
pdf(file=sprintf("%sChatNumall_%s.pdf",compare,abbr),height = 4.5, width = 4)
par(mfrow = c(1,1), xpd=TRUE)
compareInteractions(cellchat_m, show.legend = F, group = (1:2),size.text = 20)#group颜色向量 默认measure='count'
# print(compareInteractions(cellchat_m, show.legend = F, group = (1:2),size.text = 20))
compareInteractions(cellchat_m, show.legend = F, group = (1:2), measure = "weight",size.text = 20)
# print(compareInteractions(cellchat_m, show.legend = F, group = (1:2), measure = "weight",size.text = 20))
dev.off()

######不同细胞群之间的相互作用数量或强度的差异 circle
pdf(file=sprintf("%sChatNumCell_%s.pdf",compare,abbr),height = 4, width = 4)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group,sources.use=c('AP3','AP6','EC'),targets.use=c('AP3','AP6','EC'))#group为相互比较的不同的cellchat，通过names(object.list)查看用来比较的cellchat所在的位置
# print(netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group))
netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group,sources.use=c('AP3','AP6','EC'),targets.use=c('AP3','AP6','EC'))
# print(netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group))
dev.off()

####识别上调和下调的信号配体对####
pdf(file=sprintf("%sChatbubbleall_%s.pdf",compare,abbr),height = 7, width = 10)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_bubble(cellchat_m,comparison = group,sources.use=c('AP3','AP6','EC'),targets.use=c('AP3','AP6','EC'), angle.x = 90, font.size = 16) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 16))
dev.off()