library(CellChat)
library(patchwork)

load("D:/wmy/project/ALL/cellchat/cell-cell/pub.RData")
healthy<- cellchat
rm(cellchat)

load("D:/wmy/project/ALL/cellchat5_cd4cd8ya/Secreted Signaling/pub.RData")
healthy<- cellchat
rm(cellchat)
load("D:/wmy/project/ALL/cellchat5_cd4cd8ya/Secreted Signaling/no-res_post.RData")
norespost<- cellchat
rm(cellchat)
load("D:/wmy/project/ALL/cellchat5_cd4cd8ya/Secreted Signaling/no-res_pre.RData")
norespre<- cellchat
rm(cellchat)
load("D:/wmy/project/ALL/cellchat5_cd4cd8ya/Secreted Signaling/res_post.RData")
respost<- cellchat
rm(cellchat)
load("D:/wmy/project/ALL/cellchat5_cd4cd8ya/Secreted Signaling/res_pre.RData")
respre<- cellchat
rm(cellchat)

object.list <- list(healthy=healthy,respre=respre, respost = respost,norespre = norespre,norespost=norespost)
#run netAnalysis_computeCentrality
object.list<- lapply(object.list,function(x){
  x=netAnalysis_computeCentrality(x)
})
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

###取子集
names(object.list)
cellchat <- setIdent(cellchat, ident.use = "sub_clusters")
table(cellchat@idents)
cellchat
cellchat<- subsetCellChat(cellchat,idents.use =c('DPSel-7','pDC'))

######part1 从宏观角度预测细胞通讯
###比较交互总数和交互强度
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4), measure = "weight")
gg1 + gg2

###不同细胞群之间的相互作用数量或强度的差异
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T,comparison = c(1,4))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison = c(2,4))

gg1 <- netVisual_heatmap(cellchat,comparison = c(2,4))
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(2,4))
#> Do heatmap based on a merged object
gg1 + gg2


###纵览性circle plot
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, sources.use = c('DPSel-7','pDC'),targets.use = c('DPSel-7','pDC'),
                   weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

groupSize <- as.numeric(table(contrast_sec@idents))
par(mfrow = c(1,2), xpd=TRUE)
contrast_sec <- updateCellChat(contrast_sec)
table(contrast_sec@idents)
netVisual_circle(dose_sec@net$count, vertex.weight = groupSize, sources.use =c('DPSel-7','pDC') ,
                 targets.use =c('DPSel-7','pDC') ,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(dose_sec@net$weight, vertex.weight = groupSize, weight.scale = T, sources.use =c('DPSel-7','pDC') ,
                 targets.use =c('DPSel-7','pDC'),label.edge= F, title.name = "Interaction weights/strength")





########不同细胞类型之间相互作用或交互强度的差异
group.cellType <- c("DPSel-7", "pDC")
group.cellType <- factor(group.cellType, levels = c("DPSel-7", "pDC"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("weight of interactions - ", names(object.list)[i]))
}


par(mfrow = c(4,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", 
                          label.edge = T,comparison = c(2,4),title.name = 'contrast v.s. model(count)')
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", 
                          label.edge = T,comparison = c(1,2),title.name = 'contrast v.s. model(weight)')

netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", 
                          label.edge = T,comparison = c(1,3),title.name = 'contrast v.s. dose(count)')
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", 
                          label.edge = T,comparison = c(1,3),title.name = 'contrast v.s. dose(weight)')

netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", 
                          label.edge = T,comparison = c(2,3),title.name = 'model v.s. dose(count)')
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", 
                          label.edge = T,comparison = c(2,3),title.name = 'model v.s. dose(weight)')


####比较 2D 空间中的主要来源和目标
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)



pdf('D:/wmy/project/胸腺/cellchat/merge/Rplot27.pdf',18,9)
dev.off()


############################第二部分：识别保守和环境特异的信号通路##########################################
###根据信号组的功能相似性识别信号组
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional",do.parallel = F)
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)

###基于结构相似性识别信号组
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural",do.parallel = F)
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

###计算和可视化通路距离
rankSimilarity(cellchat, type = c("functional", "structural"),
               title = 'model v.s. dose',
               comparison2 = c(2,3))

#####识别并可视化保守和环境特异的信号通路
##比较每个信号通路的整体信息流
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(2, 3),)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,comparison = c(2, 3),)
gg1 + gg2

#####比较与每个细胞群相关的传出（或传入）信号
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 15)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 15)
draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 15)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 15)
draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 15)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 15)
draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))

############识别上调和下调的信号配体对
targets.use = c("pDC","Fibroblast progenitor")
table(cellchat@meta$new3)
all<- names(table(cellchat@meta$new3))
netVisual_bubble(cellchat, comparison = c(1,2,3), angle.x = 45)

all<- intersect(intersect(contrast_sec@netP$pathways,model_sec@netP$pathways),dose_sec@netP$pathways)
netVisual_bubble(cellchat, sources.use ="DPSel-7",
                 comparison = c(1,2,3), angle.x = 90 )
netVisual_bubble(cellchat, sources.use =c("pDC","Fibroblast progenitor"),targets.use = "DPSel-7",
                 comparison = c(1,2,3), angle.x = 45)

cellchat <- setIdent(cellchat, ident.use = "new3") # set "labels" as default cell identity
levels(cellchat@idents) 
###上升（增加）和下降调节（减少）信号配体受体对

all<- levels(cellchat@idents) 
tt<- c("DPbla-1","DPbla-2","DPbla-3","DPbla-4","DPbla-5","DPbla-6","DPbla-7","DPre-1", "DPre-10","DPre-2","DPre-3" ,
  "DPre-4", "DPre-5","DPre-6","DPre-7","DPre-8","DPre-9","DPSel-1","DPSel-2","DPSel-3","DPSel-4","DPSel-5",
  "DPSel-6","DPSel-7")

netVisual_bubble(cellchat, sources.use ="DPSel-7",targets.use =all,  comparison = c(2, 1), 
                        max.dataset = 2, title.name = "Increased signaling in EB", angle.x = 90, 
                        remove.isolate = T,signaling ="MHC-I")




#> Comparing communications on a merged object
netVisual_bubble(cellchat, targets.use =all,sources.use ="DPSel-7",  
                 comparison = c(2, 1), max.dataset = 2, title.name = "Increased signaling in EB", 
                 angle.x = 90, remove.isolate = T,signaling ="MHC-I")

gg1 <- netVisual_bubble(cellchat, sources.use =c("pDC","Fibroblast progenitor"),targets.use ="DPSel-7",  comparison = c(2, 3), max.dataset = 2, title.name = "Increased signaling in model", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use =c("pDC","Fibroblast progenitor"),targets.use ="DPSel-7",  comparison = c(2, 3), max.dataset = 3, title.name = "Increased signaling in dose", angle.x = 45, remove.isolate = T)

###############
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "dose"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "dose",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "dose",ligand.logFC = -0.1, receptor.logFC = -0.1)



####
raw<- cellchat
which(colnames(cellchat@net$contrast$count)=='DPSel-7')
dim(cellchat@net$contrast$count)
cellchat@net$contrast$count[1:23,1:23]<- 0
cellchat@net$model$count[1:23,1:23]<- 0
cellchat@net$dose$count[1:23,1:23]<- 0

cellchat@net$contrast$weight[1:23,1:23]<- 0
cellchat@net$model$weight[1:23,1:23]<- 0
cellchat@net$dose$weight[1:23,1:23]<- 0

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
gg1 + gg2


count_con<- cellchat@net$contrast$count
dat<- count_con[,1:3]
for (i in 1:23) {
  dat[i,1]<- count_con[24,i]+count_con[i,24]
}
dat[24,1]<- count_con[24,24]
count_mon<- cellchat@net$model$count
for (i in 1:23) {
  dat[i,2]<- count_mon[24,i]+count_mon[i,24]
}
dat[24,2]<- count_mon[24,24]
dat
count_dos<- cellchat@net$dose$count
for (i in 1:23) {
  dat[i,3]<- count_dos[24,i]+count_dos[i,24]
}
dat[24,3]<- count_dos[24,24]
dat
dat2<- data.frame(1:72,1:72,1:72)
colnames(dat2)<- c('celltype','group','count')
head(dat2)
dat2[,1]<- rep(rownames(dat),3)
dat2[,2]<- c(rep('Con',24),rep('EB',24),rep('RGLSH',24))
dat2[1:24,3]<- dat[,1]
dat2[25:48,3]<- dat[,2]
dat2[49:72,3]<- dat[,3]
library(ggplot2)
ggplot(dat2,aes(x=group,y=count,
                fill=group))+ 
  geom_bar(stat='identity')+ 
  labs(x=NULL)+             #修改坐标轴标题
  theme_bw(base_size = 18)+ #自定义主题等
  theme(axis.text = element_text(colour = 'black'),
        legend.position = 'none')


ggplot(dat2,aes(x=celltype,y=count,fill=group))+
  geom_bar(stat = 'identity', 
           #柱状图位置并排:
           position = 'dodge', #使用position=position_dodge(width=0.9),可使组内柱子间隔,自行试一下。
           width = 0.8,      #设置柱子宽度,使变量之间分开
           color='black')+        
  geom_text(aes(label=group),size=4,
            position = position_dodge(width = 0.8), #相应的注释宽度也调整
            vjust=-0.3)+    #调节注释高度    
  labs(x=NULL)+ 
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))+coord_flip()

#######
weight_con<- cellchat@net$contrast$weight
dat<- weight_con[,1:3]
for (i in 1:23) {
  dat[i,1]<- weight_con[24,i]+weight_con[i,24]
}
dat[24,1]<- weight_con[24,24]
weight_mon<- cellchat@net$model$weight
for (i in 1:23) {
  dat[i,2]<- weight_mon[24,i]+weight_mon[i,24]
}
dat[24,2]<- weight_mon[24,24]
dat
weight_dos<- cellchat@net$dose$weight
for (i in 1:23) {
  dat[i,3]<- weight_dos[24,i]+weight_dos[i,24]
}
dat[24,3]<- weight_dos[24,24]
dat
dat2<- data.frame(1:72,1:72,1:72)
colnames(dat2)<- c('celltype','group','weight')
head(dat2)
dat2[,1]<- rep(rownames(dat),3)
dat2[,2]<- c(rep('Con',24),rep('EB',24),rep('RGLSH',24))
dat2[1:24,3]<- dat[,1]
dat2[25:48,3]<- dat[,2]
dat2[49:72,3]<- dat[,3]
library(ggplot2)
ggplot(dat2,aes(x=group,y=weight,
                fill=group))+ 
  geom_bar(stat='identity')+ 
  labs(x=NULL)+             #修改坐标轴标题
  theme_bw(base_size = 18)+ #自定义主题等
  theme(axis.text = element_text(colour = 'black'),
        legend.position = 'none')


ggplot(dat2,aes(x=celltype,y=weight,fill=group))+
  geom_bar(stat = 'identity', 
           #柱状图位置并排:
           position = 'dodge', #使用position=position_dodge(width=0.9),可使组内柱子间隔,自行试一下。
           width = 0.8,      #设置柱子宽度,使变量之间分开
           color='black')+        
  geom_text(aes(label=group),size=4,
            position = position_dodge(width = 0.8), #相应的注释宽度也调整
            vjust=-0.3)+    #调节注释高度    
  labs(x=NULL)+ 
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))+coord_flip()

