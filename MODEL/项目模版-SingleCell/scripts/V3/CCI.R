####细胞通讯####
#devtools::install_github('sqjin/CellChat')
#devtools::install_github('jinworks/CellChat')最新版本
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)
# 读取数据
scRNA=readRDS('scRNA_seurat_cellType.RDS')

####1-1细胞-细胞交互####
# 拆分样本
scRNASe=subset(x = scRNA, subset = (orig.ident == "Se"))
scRNAJe=subset(x = scRNA, subset = (orig.ident == "Je"))
scRNASi=subset(x = scRNA, subset = (orig.ident == "Si"))
scRNAJi=subset(x = scRNA, subset = (orig.ident == "Ji"))
scRNAST=subset(x = scRNA, subset = (orig.ident == "ST"))
scRNAJT=subset(x = scRNA, subset = (orig.ident == "JT"))

# CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on mouse data 
# showDatabaseCategory(CellChatDB)
#预测细胞交互，并作circle图
i=1
# 样本名称
samplenames=c('Se','Si','ST','Je','Ji','JT')
# 交互方式
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

    # pdf(file=sprintf("%s%s_netVisual_circle.pdf",samplename,search),height = 15, width = 15)
    # par(mfrow = c(1,2), xpd=TRUE)
    cellchat <- updateCellChat(cellchat)
    # table(cellchat@idents)
    
    # netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
    #                 weight.scale = T, label.edge= F, title.name = "Number of interactions")

    # netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Interaction strength")
    # dev.off()
 
    # save(cellchat,file = sprintf('%s_cellchat_%s.RData',samplename,search))
    }

  i=i+1
  }

#重新绘图，进行图的调整
samplenames=c('Se','Si','ST','Je','Ji','JT')
searchs=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
for (samplename in samplenames){
  for (search in searchs){
    load(sprintf('%s_cellchat_%s.RData',samplename,search))
    # 
    # 重新绘图
    groupSize <- as.numeric(table(cellchat@idents))
    pdf(file=sprintf("%s%s_Number.pdf",samplename,search),height = 8, width = 8)
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                    weight.scale = T, label.edge= F, title.name = "Number of inferred interactions")
    dev.off()
    pdf(file=sprintf("%s%s_strength.pdf",samplename,search),height = 8, width = 8)
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Interaction strength")
    dev.off()

  }
}

# pathways 通路
load(sprintf('%s_cellchat_%s.RData','Se','Secreted Signaling'))
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "Se-net_pathway-SS.csv")

load(sprintf('%s_cellchat_%s.RData','Se','ECM-Receptor'))
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "Se-net_pathway-ECMR.csv")

load(sprintf('%s_cellchat_%s.RData','Se','Cell-Cell Contact'))
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "Se-net_pathway-CCC.csv")

#对细胞交互进行展示
# 交互方式的缩写，用于保存文件
searchsabbr=c('SS','ECM-R','CCC')
# 比较分组的样本名称，用于保存文件
comparegroup=c('SeJe','SiJi','STJT')
# 比较分组的样本在object.list所在的位置进行分组，用于后续分析
groups=list(c(1,4),c(2,5),c(3,6))
searchs=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
for (i in 1:length(searchs)){
    search=searchs[i]
    abbr=searchsabbr[i]
    # search=searchs[1]
    # abbr=searchsabbr[1]
    # 整合Se Si ST Je Ji JT Secreted Signaling cellchat
    load(sprintf('Se_cellchat_%s.RData',search))
    Se_cellchat=cellchat
    load(sprintf('Si_cellchat_%s.RData',search))
    Si_cellchat=cellchat
    load(sprintf('ST_cellchat_%s.RData',search))
    ST_cellchat=cellchat
    load(sprintf('Je_cellchat_%s.RData',search))
    Je_cellchat=cellchat
    load(sprintf('Ji_cellchat_%s.RData',search))
    Ji_cellchat=cellchat
    load(sprintf('JT_cellchat_%s.RData',search))
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
        # pdf(file=sprintf("%sChatNumall_%s.pdf",compare,abbr),height = 8, width = 8)
        # par(mfrow = c(1,1), xpd=TRUE)
        # compareInteractions(cellchat_m, show.legend = F, group = (1:6),size.text = 20)#group颜色向量 默认measure='count'
        # print(compareInteractions(cellchat_m, show.legend = F, group = (1:6),size.text = 20))
        # compareInteractions(cellchat_m, show.legend = F, group = (1:6), measure = "weight",size.text = 20)
        # print(compareInteractions(cellchat_m, show.legend = F, group = (1:6), measure = "weight",size.text = 20))
        # dev.off()
        #不同细胞群之间的相互作用数量或强度的差异 circle
        pdf(file=sprintf("%sChatNumCell_%s.pdf",compare,abbr),height = 8, width = 8)
        par(mfrow = c(1,1), xpd=TRUE)
        netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group)#group为相互比较的不同的cellchat，通过names(object.list)查看用来比较的cellchat所在的位置
        # print(netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group))
        netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group)
        # print(netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group))
        dev.off()
        #不同细胞群之间的相互作用数量或强度的差异 heatmap
        # pdf(file=sprintf("%sChatNumCellH_%s.pdf",compare,abbr),height = 8, width = 8)
        # par(mfrow = c(1,1), xpd=TRUE)
        # netVisual_heatmap(cellchat_m,comparison = group, font.size = 20, font.size.title = 20)
        # print(netVisual_heatmap(cellchat_m,comparison = group, font.size = 20, font.size.title = 20))
        # #> Do heatmap based on a merged object
        # netVisual_heatmap(cellchat_m, measure = "weight",comparison = group, font.size = 20, font.size.title = 20)
        # print(netVisual_heatmap(cellchat_m, measure = "weight",comparison = group, font.size = 20, font.size.title = 20))
        # dev.off()

        # ####识别上调和下调的信号配体对####
        # pdf(file=sprintf("%sChatbubbleall_%s.pdf",compare,abbr),height = 8, width = 15)
        # par(mfrow = c(1,1), xpd=TRUE)
        # netVisual_bubble(cellchat_m,comparison = group, angle.x = 90)
        # print(netVisual_bubble(cellchat_m,comparison = group, angle.x = 90))
        # # netVisual_bubble(cellchat_m,comparison = group, angle.x = 90)
        # # print(netVisual_bubble(cellchat_m,comparison = group, angle.x = 90))
        # dev.off()

        # pdf(file=sprintf("%sChatbubble_%s-AP.pdf",compare,abbr),height = 8, width = 8)
        # par(mfrow = c(1,1), xpd=TRUE)
        # netVisual_bubble(cellchat_m, sources.use ="AP",comparison = group, angle.x = 90)
        # print(netVisual_bubble(cellchat_m, sources.use ="AP",comparison = group, angle.x = 90))
        # netVisual_bubble(cellchat_m,targets.use = "AP",comparison = group, angle.x = 90)
        # print(netVisual_bubble(cellchat_m,targets.use = "AP",comparison = group, angle.x = 90))
        # dev.off()

        # pdf(file=sprintf("%sChatbubble_%s-myoFB.pdf",compare,abbr),height = 8, width = 8)
        # par(mfrow = c(1,1), xpd=TRUE)
        # netVisual_bubble(cellchat_m, sources.use ="myoFB",comparison = group, angle.x = 90)
        # print(netVisual_bubble(cellchat_m, sources.use ="myoFB",comparison = group, angle.x = 90))
        # netVisual_bubble(cellchat_m,targets.use = "myoFB",comparison = group, angle.x = 90)
        # print(netVisual_bubble(cellchat_m,targets.use = "myoFB",comparison = group, angle.x = 90))
        # dev.off()

        # pdf(file=sprintf("%sChatbubble_%s-APmyoFB.pdf",compare,abbr),height = 8, width = 8)
        # par(mfrow = c(1,1), xpd=TRUE)
        # netVisual_bubble(cellchat_m, sources.use ="AP",targets.use = "myoFB",comparison = group, angle.x = 90)
        # print(netVisual_bubble(cellchat_m, sources.use ="AP",targets.use = "myoFB",comparison = group, angle.x = 90))
        # netVisual_bubble(cellchat_m, sources.use ="myoFB",targets.use = "AP",comparison = group, angle.x = 90)
        # print(netVisual_bubble(cellchat_m, sources.use ="myoFB",targets.use = "AP",comparison = group, angle.x = 90))
        # dev.off()    
    }
}

####3-3细胞-细胞交互####
scRNA_SAll<- subset(scRNA,subset=(orig.ident == "Se"|orig.ident == "Si"|orig.ident == "ST"))
scRNA_JAll<- subset(scRNA,subset=(orig.ident == "Je"|orig.ident == "Ji"|orig.ident == "JT"))
i=1
samplenames=c('SAll','JAll')
searchs=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
for (seuratobj_raw in c(scRNA_SAll,scRNA_JAll)){
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
                    weight.scale = T, label.edge= F, title.name = "Number of inferred interactions")

    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,label.edge= F, title.name = "Interaction strength")
    dev.off()
 
    save(cellchat,file = sprintf('%s_cellchat_%s.RData',samplename,search))}

  i=i+1
  }

#对细胞交互进行展示
searchsabbr=c('SS','ECM-R','CCC')
comparegroup=c('SJAll')
groups=list(c(1,2))

for (i in 1:length(searchs)){
    search=searchs[i]
    abbr=searchsabbr[i]
    # search=searchs[1]
    # abbr=searchsabbr[1]
    # 整合Se Si ST Je Ji JT Secreted Signaling cellchat
    load(sprintf('SAll_cellchat_%s.RData',search))
    SAll_cellchat=cellchat
    load(sprintf('JAll_cellchat_%s.RData',search))
    JAll_cellchat=cellchat

    object.list <- list(SAll=SAll_cellchat,JAll =JAll_cellchat)
    #run netAnalysis_computeCentrality
    object.list<- lapply(object.list,function(x){
        x=netAnalysis_computeCentrality(x)})
    cellchat_m <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
    ####从宏观角度预测细胞通讯####
    #比较交互总数和交互强度

# for (num in 1:length(searchsabbr)){
    compare=comparegroup[1]
    group=groups[[1]]
    # compare=comparegroup[3]
    # group=groups[[3]]
    pdf(file=sprintf("%sChatNumall_%s.pdf",compare,abbr),height = 8, width = 8)
    par(mfrow = c(1,1), xpd=TRUE)
    compareInteractions(cellchat_m, show.legend = F, group = (1:2),size.text = 20)#group颜色向量
    print(compareInteractions(cellchat_m, show.legend = F, group = (1:2),size.text = 20))
    compareInteractions(cellchat_m, show.legend = F, group = (1:2), measure = "weight",size.text = 20)
    print(compareInteractions(cellchat_m, show.legend = F, group = (1:2), measure = "weight",size.text = 20))
    dev.off()
    #不同细胞群之间的相互作用数量或强度的差异 circle
    pdf(file=sprintf("%sChatNumCell_%s.pdf",compare,abbr),height = 8, width = 8)
    par(mfrow = c(1,1), xpd=TRUE)
    netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group)#group为相互比较的不同的cellchat，通过names(object.list)查看用来比较的cellchat所在的位置
    # print(netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group))
    netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group)
    # print(netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group))
    dev.off()
    #不同细胞群之间的相互作用数量或强度的差异 heatmap
    pdf(file=sprintf("%sChatNumCellH_%s.pdf",compare,abbr),height = 8, width = 8)
    par(mfrow = c(1,1), xpd=TRUE)
    netVisual_heatmap(cellchat_m,comparison = group, font.size = 20, font.size.title = 20)
    print(netVisual_heatmap(cellchat_m,comparison = group, font.size = 20, font.size.title = 20))
    #> Do heatmap based on a merged object
    netVisual_heatmap(cellchat_m, measure = "weight",comparison = group, font.size = 20, font.size.title = 20)
    print(netVisual_heatmap(cellchat_m, measure = "weight",comparison = group, font.size = 20, font.size.title = 20))
    dev.off()

    ####识别上调和下调的信号配体对####
    pdf(file=sprintf("%sChatbubble_%s.pdf",compare,abbr),height = 8, width = 15)
    par(mfrow = c(1,1), xpd=TRUE)
    netVisual_bubble(cellchat_m,comparison = group, angle.x = 90)
    print(netVisual_bubble(cellchat_m,comparison = group, angle.x = 90))
    netVisual_bubble(cellchat_m,comparison = group, angle.x = 90)
    print(netVisual_bubble(cellchat_m,comparison = group, angle.x = 90))
    dev.off()

    pdf(file=sprintf("%sChatbubble_%s-AP.pdf",compare,abbr),height = 8, width = 8)
    par(mfrow = c(1,1), xpd=TRUE)
    netVisual_bubble(cellchat_m, sources.use ="AP",comparison = group, angle.x = 90)
    print(netVisual_bubble(cellchat_m, sources.use ="AP",comparison = group, angle.x = 90))
    netVisual_bubble(cellchat_m,targets.use = "AP",comparison = group, angle.x = 90)
    print(netVisual_bubble(cellchat_m,targets.use = "AP",comparison = group, angle.x = 90))
    dev.off()
# }
}


# 单个样本配受体-特定细胞可视化
searchsabbr=c('SS','ECM-R','CCC')
search='Cell-Cell Contact'
for (search in searchs){
  # searchs=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
  load(sprintf('Se_cellchat_%s.RData',search))
  levels(cellchat@idents)
  # show all the significant interactions (L-R pairs)
  #需要指定受体细胞和配体细胞
  # p = netVisual_bubble(cellchat, remove.isolate = FALSE, angle.x = 90)
  # ggsave("Se_bubbleSingle-%s.pdf", p, width = 12, height = 8) 
  p = netVisual_bubble(cellchat, sources.use = c('AP','myoFB','EC'), 
                      targets.use = c('AP','myoFB','EC'), remove.isolate = FALSE,angle.x = 45)+theme(axis.text.x=element_text(vjust=1,size=16,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=16,face = "bold"))
    
  ggsave(sprintf("Se_bubbleSingleAPmyoFBEC-%s.pdf",search), p, width = 6, height = 6)
  load(sprintf('Je_cellchat_%s.RData',search))
  levels(cellchat@idents)
  # show all the significant interactions (L-R pairs)
  #需要指定受体细胞和配体细胞
  # p = netVisual_bubble(cellchat, remove.isolate = FALSE, angle.x = 90)
  # ggsave("Je_bubbleSingle-%s.pdf", p, width = 12, height = 8)
  p = netVisual_bubble(cellchat, sources.use = c('AP','myoFB','EC'), 
                      targets.use = c('AP','myoFB','EC'), remove.isolate = FALSE,angle.x = 45)+theme(axis.text.x=element_text(vjust=1,size=16,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=16,face = "bold"))
    
  ggsave(sprintf("Je_bubbleSingleAPmyoFBEC-%s.pdf",search), p, width = 6, height = 6)
}


####3组3模式相互作用数量及强度展示
searchs=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
searchsabbr=c('SS','ECM-R','CCC')
comparegroup=c('SeJe','SiJi','STJT')
# 比较分组的样本在object.list所在的位置进行分组，用于后续分析
# groups=list(c(1,4),c(2,5),c(3,6))

for (num in 1:length(searchs)){
  search=searchs[num]
  abbr=searchsabbr[num]
  # 整合Se Si ST Je Ji JT Secreted Signaling cellchat
  load(sprintf('Se_cellchat_%s.RData',search))
  Se_cellchat=cellchat
  load(sprintf('Si_cellchat_%s.RData',search))
  Si_cellchat=cellchat
  load(sprintf('ST_cellchat_%s.RData',search))
  ST_cellchat=cellchat
  load(sprintf('Je_cellchat_%s.RData',search))
  Je_cellchat=cellchat
  load(sprintf('Ji_cellchat_%s.RData',search))
  Ji_cellchat=cellchat
  load(sprintf('JT_cellchat_%s.RData',search))
  JT_cellchat=cellchat

  # object.list <- list(Se=Se_cellchat,Si =Si_cellchat,ST=ST_cellchat,Je=Je_cellchat,Ji=Ji_cellchat,JT=JT_cellchat)
  i=3
  if(i==1) {
    object.list <- list(Se=Se_cellchat,Je=Je_cellchat)#i=1
  } else if (i==2) {
    object.list <- list(Si=Si_cellchat,Ji=Ji_cellchat)#i=2
  } else if (i==3){
    object.list <- list(ST=ST_cellchat,JT=JT_cellchat)#i=3
  }
  # object.list <- list(Se=Se_cellchat,Je=Je_cellchat)#i=1
  # object.list <- list(Si=Si_cellchat,Ji=Ji_cellchat)#i=2
  # object.list <- list(ST=ST_cellchat,JT=JT_cellchat)#i=3
  #run netAnalysis_computeCentrality
  object.list<- lapply(object.list,function(x){
      x=netAnalysis_computeCentrality(x)})
  cellchat_m <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
  ####从宏观角度预测细胞通讯####
  #比较交互总数和交互强度

  compare=comparegroup[i]
  group=groups[[i]]
  pdf(file=sprintf("%sChatNumall_%s.pdf",compare,abbr),height = 4, width = 4)
  par(mfrow = c(1,2), xpd=TRUE)
  compareInteractions(cellchat_m, show.legend = F, group = (1:2),size.text = 18)
  print(compareInteractions(cellchat_m, show.legend = F, group = (1:2),size.text = 18))
  compareInteractions(cellchat_m, show.legend = F, group = (1:2), measure = "weight",size.text = 18)
  print(compareInteractions(cellchat_m, show.legend = F, group = (1:2), measure = "weight",size.text = 18))
  dev.off()
}

# Se-Je-3种模式下细胞间的 AP差异数量及强度展示
searchs=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
searchsabbr=c('SS','ECM-R','CCC')
# comparegroup=c('SeJe','SiJi','STJT')
for (num in 1:length(searchs)){
  search=searchs[num]
  abbr=searchsabbr[num]
  # 整合Se Si ST Je Ji JT Secreted Signaling cellchat
  load(sprintf('Se_cellchat_%s.RData',search))
  Se_cellchat=cellchat
  # load(sprintf('Si_cellchat_%s.RData',search))
  # Si_cellchat=cellchat
  # load(sprintf('ST_cellchat_%s.RData',search))
  # ST_cellchat=cellchat
  load(sprintf('Je_cellchat_%s.RData',search))
  Je_cellchat=cellchat
  load(sprintf('Ji_cellchat_%s.RData',search))
  # Ji_cellchat=cellchat
  # load(sprintf('JT_cellchat_%s.RData',search))
  # JT_cellchat=cellchat
  i = 1#通过rownames(mat_se_count)查看AP位置index
  
  mat_se_count <- Se_cellchat@net$count
  mat_se_count2 <- matrix(0, nrow = nrow(mat_se_count), ncol = ncol(mat_se_count), dimnames = dimnames(mat_se_count))
  mat_se_count2[i, ] <- mat_se_count[i, ]
  mat_se_weight <- Se_cellchat@net$weight
  mat_se_weight2 <- matrix(0, nrow = nrow(mat_se_weight), ncol = ncol(mat_se_weight), dimnames = dimnames(mat_se_weight))
  mat_se_weight2[i, ] <- mat_se_weight[i, ]
  groupSize <- as.numeric(table(Se_cellchat@idents))
  pdf(file=sprintf("Se-NumStrength_%s.pdf",abbr),height = 4, width = 4)
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle(mat_se_count2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                    arrow.size = 0.1, edge.weight.max = max(mat_se_count), title.name = 'Number of inferred interactions')
  netVisual_circle(mat_se_weight2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                    arrow.size = 0.1, edge.weight.max = max(mat_se_weight), title.name = 'Interaction strength')  
  dev.off()
  pdf(file=sprintf("Se-NumStrengthdiff_%s.pdf",abbr),height = 4, width = 4)
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group)
  dev.off()

  mat_je_count <- Je_cellchat@net$count
  mat_je_count2 <- matrix(0, nrow = nrow(mat_je_count), ncol = ncol(mat_je_count), dimnames = dimnames(mat_je_count))
  mat_je_count2[i, ] <- mat_je_count[i, ]
  mat_je_weight <- Je_cellchat@net$weight
  mat_je_weight2 <- matrix(0, nrow = nrow(mat_je_weight), ncol = ncol(mat_je_weight), dimnames = dimnames(mat_je_weight))
  mat_je_weight2[i, ] <- mat_je_weight[i, ]
  groupSize <- as.numeric(table(Je_cellchat@idents))
  pdf(file=sprintf("Je-NumStrength_%s.pdf",abbr),height = 4, width = 4)
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle(mat_je_count2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                    arrow.size = 0.1, edge.weight.max = max(mat_je_count), title.name = 'Number of inferred interactions')
  netVisual_circle(mat_je_weight2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                    arrow.size = 0.1, edge.weight.max = max(mat_je_weight), title.name = 'Interaction strength')  
                    
  dev.off()
}


# 1.加载Seurat数据+设置工作目录----#####

Se_gene=unique(c("Sema3c","Nampt","Igf1","Gas6","Fgr10","Adipoq","Nrp1","Plxna2","Insr","Mertk","Fgfr2","Fgfr1","Adipor2","Sema3c","Insr","Mertk","Fgfr2","Adipor2","Lamc1","Lamb1","Lama4","Dag1","Itga7","Itga6","Itgb1","Lama1","Dag1","Itga6","Itgb1","Ptprm","Negr1","Jag1","Ptprm","Negr1","Notch2","Fgfr1","Nrp1","Plxna2","Sema4a","Ptprm","Ncam1","Mpzl1","Cdh1","Cadm1","Ptprm","Jag1","Cdh1","Ncam1","Mpzl1","Cdh1","Cadm1"))

exprSet <- scRNASe@assays$RNA@data
exprSet<-as.data.frame(t(exprSet))#转置

# 2. 计算某个基因和其它基因的相关性（以S100A8为例）-----#####
y <- as.numeric(exprSet[,"Nucb2"])

col_pos <- match(Se_gene, colnames(exprSet))
col_pos=col_pos[!is.na(col_pos)]
df_existing <- exprSet[,col_pos]
colnames<-colnames(df_existing)
cor_data_df<- data.frame(colnames)
for(i in 1:length(colnames)){
  test_p<-cor.test(as.numeric(df_existing[,i]),y,method = "pearson")
  test_s<-cor.test(as.numeric(df_existing[,i]),y,method = "spearman")
  test_k<-cor.test(as.numeric(df_existing[,i]),y,method = "kendall")
  cor_data_df[i,2]<- test_s$estimate
  cor_data_df[i,3]<- test_s$p.value 
  cor_data_df[i,4]<- test_p$estimate
  cor_data_df[i,5]<- test_p$p.value
  cor_data_df[i,6]<- test_k$estimate
  cor_data_df[i,7]<- test_k$p.value 
}
names(cor_data_df)<-c("symbol","correlation-s","pvalue-s","correlation-p","pvalue-p","correlation-k","pvalue-k")
# cor_data_df %>% head()
write.csv(file='Se.Nucb2-interaction.cor.csv',cor_data_df)

# 3. 筛选有意义的正相关和负相关的基因-----####
library(dplyr)
library(tidyr)
cor_data_sig_pos <- cor_data_df %>%
  dplyr::filter(pvalue <0.01)%>%dplyr::filter(correlation >0)%>%
  dplyr::arrange(desc(correlation))
 
cor_data_sig_neg <- cor_data_df %>%
  dplyr::filter(pvalue <0.01)%>%dplyr::filter(correlation <0)%>%
  dplyr::arrange(desc(abs(correlation)))

# ----------------------------------------------------------------
Je_gene=unique(c("Sema3c","Retn","Rarres2","Nrg4","Nampt","Fgf2","Fgf10","Angptl2","Adipoq","Nrp1","Plxna2","Cmklr1","Erbb4","Insr","Fgfr2","Fgfr1","Adipor2","lamc1","lamb1","lama4","Col4a2","Dag1","Itga7","Itgb1","Itga1","Ptprm","Negr1","Jag1","Cd46","Ptprm","Negr1","Jag1","Notch2"))
exprSet <- scRNAJe@assays$RNA@data
exprSet<-as.data.frame(t(exprSet))#转置

# 2. 计算某个基因和其它基因的相关性（以S100A8为例）-----#####
y <- as.numeric(exprSet[,"Nucb2"])

col_pos <- match(Je_gene, colnames(exprSet))
col_pos=col_pos[!is.na(col_pos)]
df_existing <- exprSet[,col_pos]
colnames<-colnames(df_existing)
cor_data_df<- data.frame(colnames)
for(i in 1:length(colnames)){
  test_p<-cor.test(as.numeric(df_existing[,i]),y,method = "pearson")
  test_s<-cor.test(as.numeric(df_existing[,i]),y,method = "spearman")
  test_k<-cor.test(as.numeric(df_existing[,i]),y,method = "kendall")
  cor_data_df[i,2]<- test_s$estimate
  cor_data_df[i,3]<- test_s$p.value 
  cor_data_df[i,4]<- test_p$estimate
  cor_data_df[i,5]<- test_p$p.value
  cor_data_df[i,6]<- test_k$estimate
  cor_data_df[i,7]<- test_k$p.value 
}
names(cor_data_df)<-c("symbol","correlation-s","pvalue-s","correlation-p","pvalue-p","correlation-k","pvalue-k")
# cor_data_df %>% head()
write.csv(file='Je.Nucb2-interaction.cor.csv',cor_data_df)
# 3. 筛选有意义的正相关和负相关的基因-----####
library(dplyr)
library(tidyr)
cor_data_sig_pos <- cor_data_df %>%
  dplyr::filter(pvalue <0.01)%>%dplyr::filter(correlation >0)%>%
  dplyr::arrange(desc(correlation))
 
cor_data_sig_neg <- cor_data_df %>%
  dplyr::filter(pvalue <0.01)%>%dplyr::filter(correlation <0)%>%
  dplyr::arrange(desc(abs(correlation)))


# 4. 随机选取正相关和负相关基因，分别作图验证----######
#1）S100A9正相关####
# install.packages("ggstatsplot") ##linux安装过程中报错，需先安装mpfr工具
# $conda install anaconda::mpfr
library(ggstatsplot)
png(paste0("相关性densigram-正相关.png"), width = 6, height = 6, res = 400, units = "in")
ggscatterstats(data = exprSet ,
               y = S100A8,
               x = S100A9,
               centrality.para = "mean",               
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               marginal.type = "densigram", # #类型可以换成density,boxplot,violin,densigram                
               title = "Relationship between S100A8 and S100A9")
dev.off()
#2）S100A9负相关####
library(ggstatsplot)
png(paste0("相关性densigram-负相关.png"), width = 6, height = 6, res = 400, units = "in")
ggscatterstats(data = exprSet,
               y = MALAT1,
               x = S100A9,
               centrality.para = "mean",               
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               marginal.type = "densigram", # #类型可以换成density,boxplot,violin,densigram                
               title = "Relationship between MALAT1 and S100A9")
dev.off()
 
# 这两个图的更简单画法（Seurat的FeatureScatter函数)
library(Seurat)
library(SeuratObject)
pdf(file=sprintf("Se-Nucb2.cor.pdf",abbr),height = 4, width = 6)
par(mfrow = c(2,3), xpd=TRUE)

FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Sema3c")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Ncam1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Lama1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Fgfr2")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Igf1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Nrp1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))

dev.off()

Idents(scRNASe) <- 'cellType'
#去除Nucb2表达为0的细胞
filtered_cells <- WhichCells(scRNASe, expression = Nucb2 > 0)
seurat_n <- subset(scRNASe, cells = filtered_cells)
pdf(file=sprintf("Se-Nucb2-filt.cor.pdf"),height = 4, width = 6)
par(mfrow = c(2,3), xpd=TRUE)

filtered_cells <- WhichCells(seurat_n, expression = Sema3c > 0)
seurat_ns <- subset(seurat_n, cells = filtered_cells)
FeatureScatter(seurat_ns, feature1 = "Nucb2", feature2 = "Sema3c")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
filtered_cells <- WhichCells(seurat_n, expression = Ncam1 > 0)
seurat_nc <- subset(seurat_n, cells = filtered_cells)
FeatureScatter(seurat_nc, feature1 = "Nucb2", feature2 = "Ncam1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
filtered_cells <- WhichCells(seurat_n, expression = Lama1 > 0)
seurat_nl <- subset(seurat_n, cells = filtered_cells)
FeatureScatter(seurat_nl, feature1 = "Nucb2", feature2 = "Lama1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
filtered_cells <- WhichCells(seurat_n, expression = Fgfr2 > 0)
seurat_nf <- subset(seurat_n, cells = filtered_cells)
FeatureScatter(seurat_nf, feature1 = "Nucb2", feature2 = "Fgfr2")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
filtered_cells <- WhichCells(seurat_n, expression = Igf1 > 0)
seurat_ni <- subset(seurat_n, cells = filtered_cells)
FeatureScatter(seurat_ni,feature1 = "Nucb2", feature2 = "Igf1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
filtered_cells <- WhichCells(seurat_n, expression = Nrp1 > 0)
seurat_nr <- subset(seurat_n, cells = filtered_cells)
FeatureScatter(seurat_nr, feature1 = "Nucb2", feature2 = "Nrp1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))

dev.off()



scRNASeAP=subset(scRNASe,subset=(cellType =="AP"))
scRNASemyoFB=subset(scRNASe,subset=(cellType =="myoFB"))

pdf(file=sprintf("Se-AP-Nucb2.cor.pdf",abbr),height = 4, width = 6)
par(mfrow = c(2,3), xpd=TRUE)

FeatureScatter(scRNASeAP, feature1 = "Nucb2", feature2 = "Sema3c")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASeAP, feature1 = "Nucb2", feature2 = "Ncam1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASeAP, feature1 = "Nucb2", feature2 = "Lama1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASeAP, feature1 = "Nucb2", feature2 = "Fgfr2")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASeAP, feature1 = "Nucb2", feature2 = "Igf1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASeAP, feature1 = "Nucb2", feature2 = "Nrp1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))

dev.off()


myoFBFeature=unique(c("Sema3c","Vegfa","Insr","Mertk","Fgfr2","Adipor2","Igf1r","Lama1","Col6a4","Dag1","Itga6","Itgb1","Sema4a","Ptprm","Ncam1","Mpzl1","Cdh1","Cadm1","App","Ptprm","Jag1","Cdh1","Ncam1","Mpzl1","Cdh1","Cadm1","Notch2","Notch1"))
exprSet <- scRNASemyoFB@assays$RNA@data
exprSet<-as.data.frame(t(exprSet))#转置

# 2. 计算某个基因和其它基因的相关性（以S100A8为例）-----#####
y <- as.numeric(exprSet[,"Nucb2"])

col_pos <- match(myoFBFeature, colnames(exprSet))
col_pos=col_pos[!is.na(col_pos)]
df_existing <- exprSet[,col_pos]
colnames<-colnames(df_existing)
cor_data_df<- data.frame(colnames)
for(i in 1:length(colnames)){
  test_p<-cor.test(as.numeric(df_existing[,i]),y,method = "pearson")
  test_s<-cor.test(as.numeric(df_existing[,i]),y,method = "spearman")
  test_k<-cor.test(as.numeric(df_existing[,i]),y,method = "kendall")
  cor_data_df[i,2]<- test_s$estimate
  cor_data_df[i,3]<- test_s$p.value 
  cor_data_df[i,4]<- test_p$estimate
  cor_data_df[i,5]<- test_p$p.value
  cor_data_df[i,6]<- test_k$estimate
  cor_data_df[i,7]<- test_k$p.value 
}
names(cor_data_df)<-c("symbol","correlation-s","pvalue-s","correlation-p","pvalue-p","correlation-k","pvalue-k")
# cor_data_df %>% head()
write.csv(file='SemyoFBFeature.Nucb2-interaction.cor.csv',cor_data_df)

pdf(file=sprintf("Se-myoFB-Nucb2.cor.pdf",height = 4, width = 6))
par(mfrow = c(2,3), xpd=TRUE)

FeatureScatter(scRNASemyoFB, feature1 = "Nucb2", feature2 = "Sema3c")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASemyoFB, feature1 = "Nucb2", feature2 = "Ptprm")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASemyoFB, feature1 = "Nucb2", feature2 = "Lama1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASemyoFB, feature1 = "Nucb2", feature2 = "Igf1r")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))

dev.off()

pdf(file=sprintf("Se-myoFB-Nucb2.cor-Igf.pdf"),height = 4, width = 6)
par(mfrow = c(2,3), xpd=TRUE)
FeatureScatter(scRNASemyoFB, feature1 = "Nucb2", feature2 = "Igf1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASemyoFB, feature1 = "Nucb2", feature2 = "Igf1r")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
dev.off()

pdf(file=sprintf("Se-EC-Nucb2.cor-VEGF.pdf"),height = 4, width = 6)
# Vegfa，Vegfb，Vegfr1，Vegfr2，Vegfr1r2
par(mfrow = c(2,3), xpd=TRUE)
FeatureScatter(scRNASeEC, feature1 = "Nucb2", feature2 = "Vegfa")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASeEC, feature1 = "Nucb2", feature2 = "Vegfb")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASeEC, feature1 = "Nucb2", feature2 = "Vegfr1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASeEC, feature1 = "Nucb2", feature2 = "Vegfr2")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASeEC, feature1 = "Nucb2", feature2 = "Vegfr1r2")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
dev.off()

# AP为中心的差异配受体数量和强度
searchsabbr=c('SS','ECM-R','CCC')
# 比较分组的样本名称，用于保存文件
comparegroup=c('SeJe','SiJi','STJT')
# 比较分组的样本在object.list所在的位置进行分组，用于后续分析
groups=list(c(4,1),c(5,2),c(6,3))
searchs=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
for (i in 1:length(searchs)){
    search=searchs[i]
    abbr=searchsabbr[i]
    # search=searchs[1]
    # abbr=searchsabbr[1]
    # 整合Se Si ST Je Ji JT Secreted Signaling cellchat
    load(sprintf('Se_cellchat_%s.RData',search))
    Se_cellchat=cellchat
    load(sprintf('Si_cellchat_%s.RData',search))
    Si_cellchat=cellchat
    load(sprintf('ST_cellchat_%s.RData',search))
    ST_cellchat=cellchat
    load(sprintf('Je_cellchat_%s.RData',search))
    Je_cellchat=cellchat
    load(sprintf('Ji_cellchat_%s.RData',search))
    Ji_cellchat=cellchat
    load(sprintf('JT_cellchat_%s.RData',search))
    JT_cellchat=cellchat

    object.list <- list(Se=Se_cellchat,Si =Si_cellchat,ST=ST_cellchat,Je=Je_cellchat,Ji=Ji_cellchat,JT=JT_cellchat)
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

    }
}

Idents(scRNASe) <- 'cellType'
#去除Nucb2表达为0的细胞
filtered_cells <- WhichCells(scRNASe, expression = Nucb2 > 0)
seurat_n <- subset(scRNASe, cells = filtered_cells)
Se_gene=unique(c("Sema3c","Nampt","Igf1","Gas6","Fgr10","Adipoq","Nrp1","Plxna2","Insr","Mertk","Fgfr2","Fgfr1","Adipor2","Sema3c","Insr","Mertk","Fgfr2","Adipor2","Lamc1","Lamb1","Lama4","Dag1","Itga7","Itga6","Itgb1","Lama1","Dag1","Itga6","Itgb1","Ptprm","Negr1","Jag1","Ptprm","Negr1","Notch2","Fgfr1","Nrp1","Plxna2","Sema4a","Ptprm","Ncam1","Mpzl1","Cdh1","Cadm1","Ptprm","Jag1","Cdh1","Ncam1","Mpzl1","Cdh1","Cadm1","Vegfa","Igf1r","Gas6","Vegfr1","Nrp1","Plxna2","Nrp2","Plxna4","Insr","Igf1r","Fgfr2","Adipor2","Col4a2","Col4a1","Col6a4","Lamc1","Lamb1","Lama4","Hspg2","Col4a2","Col4a1","Col1a2","Itga6","Itgb1","Itga1"))

# "Vegfa","Igf1r","Gas6","Vegfr1","Nrp1","Plxna2","Nrp2","Plxna4","Insr","Igf1r","Fgfr2","Adipor2","Col4a2","Col4a1","Col6a4","Lamc1","Lamb1","Lama4","Hspg2","Col4a2","Col4a1","Col1a2","Itga6","Itgb1","Itga1"

exprSet <- seurat_n@assays$RNA@data
exprSet<-as.data.frame(t(exprSet))#转置

# 2. 计算某个基因和其它基因的相关性#####
y <- as.numeric(exprSet[,"Nucb2"])

col_pos <- match(Se_gene, colnames(exprSet))
col_pos=col_pos[!is.na(col_pos)]
df_existing <- exprSet[,col_pos]
colnames<-colnames(df_existing)
cor_data_df<- data.frame(colnames)
for(i in 1:length(colnames)){
  test_p<-cor.test(as.numeric(df_existing[,i]),y,method = "pearson")
  test_s<-cor.test(as.numeric(df_existing[,i]),y,method = "spearman")
  test_k<-cor.test(as.numeric(df_existing[,i]),y,method = "kendall")
  cor_data_df[i,2]<- test_s$estimate
  cor_data_df[i,3]<- test_s$p.value 
  cor_data_df[i,4]<- test_p$estimate
  cor_data_df[i,5]<- test_p$p.value
  cor_data_df[i,6]<- test_k$estimate
  cor_data_df[i,7]<- test_k$p.value 
}
names(cor_data_df)<-c("symbol","correlation-s","pvalue-s","correlation-p","pvalue-p","correlation-k","pvalue-k")
# cor_data_df %>% head()
write.csv(file='Se.Nucb2-interaction.cor2.csv',cor_data_df)

Idents(scRNAJe) <- 'cellType'
#去除Nucb2表达为0的细胞
filtered_cells <- WhichCells(scRNAJe, expression = Nucb2 > 0)
seurat_n <- subset(scRNAJe, cells = filtered_cells)
Je_gene=unique(c("Sema3c","Retn","Rarres2","Nrg4","Nampt","Fgf2","Fgf10","Angptl2","Adipoq","Nrp1","Plxna2","Cmklr1","Erbb4","Insr","Fgfr2","Fgfr1","Adipor2","lamc1","lamb1","lama4","Col4a2","Dag1","Itga7","Itgb1","Itga1","Ptprm","Negr1","Jag1","Cd46","Ptprm","Negr1","Jag1","Notch2","Vegfb","Vegfa","Igf1","Angptl4","Angpt1","Tlr4","Bmpr1a","Bmpr2","Acvr2a","Acvr1","Adipoq","Bmp6","Retn","Vegfr1","Vegfr2","Vegfr1r2","Nrp1","Plxna2","Nrp2","Plxna4","Cdh5","Tek","Adipor2","Col4a1","Col6a3","Lamc1","Lama4","Hspg2","Col4a2","Col1a2","Itga6","Itgb1","Itga1","App","Plxna2","Itgav","Itgb1","Sema6a","Ptprm","Pecam1","Jam2","Esam","Cdh5","APP","Ptprm","Cd74","Plxna4","Plxna2","Ptprm","Pecam1","Jam2","Esam","Cdh5","Cd74"))

# "Vegfb","Vegfa","Igf1","Angptl4","Angpt1","Tlr4","Bmpr1a","Bmpr2","Acvr2a","Acvr1","Adipoq","Bmp6","Retn","Vegfr1","Vegfr2","Vegfr1r2","Nrp1","Plxna2","Nrp2","Plxna4","Cdh5","Tek","Adipor2","Col4a1","Col6a3","Lamc1","Lama4","Hspg2","Col4a2","Col1a2","Itga6","Itgb1","Itga1","App","Plxna2","Itgav","Itgb1","Sema6a","Ptprm","Pecam1","Jam2","Esam","Cdh5","APP","Ptprm","Cd74","Plxna4","Plxna2","Ptprm","Pecam1","Jam2","Esam","Cdh5","Cd74"

exprSet <- seurat_n@assays$RNA@data
exprSet<-as.data.frame(t(exprSet))#转置

# 2. 计算某个基因和其它基因的相关性
y <- as.numeric(exprSet[,"Nucb2"])

col_pos <- match(Je_gene, colnames(exprSet))
col_pos=col_pos[!is.na(col_pos)]
df_existing <- exprSet[,col_pos]
colnames<-colnames(df_existing)
cor_data_df<- data.frame(colnames)
for(i in 1:length(colnames)){
  test_p<-cor.test(as.numeric(df_existing[,i]),y,method = "pearson")
  test_s<-cor.test(as.numeric(df_existing[,i]),y,method = "spearman")
  test_k<-cor.test(as.numeric(df_existing[,i]),y,method = "kendall")
  cor_data_df[i,2]<- test_s$estimate
  cor_data_df[i,3]<- test_s$p.value 
  cor_data_df[i,4]<- test_p$estimate
  cor_data_df[i,5]<- test_p$p.value
  cor_data_df[i,6]<- test_k$estimate
  cor_data_df[i,7]<- test_k$p.value 
}
names(cor_data_df)<-c("symbol","correlation-s","pvalue-s","correlation-p","pvalue-p","correlation-k","pvalue-k")
# cor_data_df %>% head()
write.csv(file='Je.Nucb2-interaction.cor2.csv',cor_data_df)

library(Seurat)
library(SeuratObject)
pdf(file=sprintf("Se-Nucb2.cor2.pdf"),height = 4, width = 6)
par(mfrow = c(2,3), xpd=TRUE)
FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Vegfa")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Notch2")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Insr")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Col6a4")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Itga6")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Cdh1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Mpzl1")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
FeatureScatter(scRNASe, feature1 = "Nucb2", feature2 = "Igf1r")+ theme(text = element_text(size = 20))+theme(axis.text.x=element_text(vjust=1,size=10,face = "bold"))+theme(axis.text.y=element_text(vjust=1,size=10,face = "bold"))
dev.off()


Idents(scRNASe) <- 'cellType'
#去除Nucb2表达为0的细胞
# filtered_cells <- WhichCells(scRNASe, expression = Nucb2 > 0)
# seurat_n <- subset(scRNASe, cells = filtered_cells)
Se_gene=unique(c("Sema3c","Nampt","Igf1","Gas6","Fgr10","Adipoq","Nrp1","Plxna2","Insr","Mertk","Fgfr2","Fgfr1","Adipor2","Sema3c","Insr","Mertk","Fgfr2","Adipor2","Lamc1","Lamb1","Lama4","Dag1","Itga7","Itga6","Itgb1","Lama1","Dag1","Itga6","Itgb1","Ptprm","Negr1","Jag1","Ptprm","Negr1","Notch2","Fgfr1","Nrp1","Plxna2","Sema4a","Ptprm","Ncam1","Mpzl1","Cdh1","Cadm1","Ptprm","Jag1","Cdh1","Ncam1","Mpzl1","Cdh1","Cadm1","Vegfa","Igf1r","Gas6","Vegfr1","Nrp1","Plxna2","Nrp2","Plxna4","Insr","Igf1r","Fgfr2","Adipor2","Col4a2","Col4a1","Col6a4","Lamc1","Lamb1","Lama4","Hspg2","Col4a2","Col4a1","Col1a2","Itga6","Itgb1","Itga1"))

# "Vegfa","Igf1r","Gas6","Vegfr1","Nrp1","Plxna2","Nrp2","Plxna4","Insr","Igf1r","Fgfr2","Adipor2","Col4a2","Col4a1","Col6a4","Lamc1","Lamb1","Lama4","Hspg2","Col4a2","Col4a1","Col1a2","Itga6","Itgb1","Itga1"

exprSet <- scRNASe@assays$RNA@data
exprSet<-as.data.frame(t(exprSet))#转置

# 2. 计算某个基因和其它基因的相关性#####
y <- as.numeric(exprSet[,"Nucb2"])

col_pos <- match(Se_gene, colnames(exprSet))
col_pos=col_pos[!is.na(col_pos)]
df_existing <- exprSet[,col_pos]
colnames<-colnames(df_existing)
cor_data_df<- data.frame(colnames)
for(i in 1:length(colnames)){
  test_p<-cor.test(as.numeric(df_existing[,i]),y,method = "pearson")
  test_s<-cor.test(as.numeric(df_existing[,i]),y,method = "spearman")
  test_k<-cor.test(as.numeric(df_existing[,i]),y,method = "kendall")
  cor_data_df[i,2]<- test_s$estimate
  cor_data_df[i,3]<- test_s$p.value 
  cor_data_df[i,4]<- test_p$estimate
  cor_data_df[i,5]<- test_p$p.value
  cor_data_df[i,6]<- test_k$estimate
  cor_data_df[i,7]<- test_k$p.value 
}
names(cor_data_df)<-c("symbol","correlation-s","pvalue-s","correlation-p","pvalue-p","correlation-k","pvalue-k")
# cor_data_df %>% head()
write.csv(file='Se.Nucb2-interaction.cor3.csv',cor_data_df)

Idents(scRNAJe) <- 'cellType'
#去除Nucb2表达为0的细胞
# filtered_cells <- WhichCells(scRNAJe, expression = Nucb2 > 0)
# seurat_n <- subset(scRNAJe, cells = filtered_cells)
Je_gene=unique(c("Sema3c","Retn","Rarres2","Nrg4","Nampt","Fgf2","Fgf10","Angptl2","Adipoq","Nrp1","Plxna2","Cmklr1","Erbb4","Insr","Fgfr2","Fgfr1","Adipor2","lamc1","lamb1","lama4","Col4a2","Dag1","Itga7","Itgb1","Itga1","Ptprm","Negr1","Jag1","Cd46","Ptprm","Negr1","Jag1","Notch2","Vegfb","Vegfa","Igf1","Angptl4","Angpt1","Tlr4","Bmpr1a","Bmpr2","Acvr2a","Acvr1","Adipoq","Bmp6","Retn","Vegfr1","Vegfr2","Vegfr1r2","Nrp1","Plxna2","Nrp2","Plxna4","Cdh5","Tek","Adipor2","Col4a1","Col6a3","Lamc1","Lama4","Hspg2","Col4a2","Col1a2","Itga6","Itgb1","Itga1","App","Plxna2","Itgav","Itgb1","Sema6a","Ptprm","Pecam1","Jam2","Esam","Cdh5","APP","Ptprm","Cd74","Plxna4","Plxna2","Ptprm","Pecam1","Jam2","Esam","Cdh5","Cd74"))

# "Vegfb","Vegfa","Igf1","Angptl4","Angpt1","Tlr4","Bmpr1a","Bmpr2","Acvr2a","Acvr1","Adipoq","Bmp6","Retn","Vegfr1","Vegfr2","Vegfr1r2","Nrp1","Plxna2","Nrp2","Plxna4","Cdh5","Tek","Adipor2","Col4a1","Col6a3","Lamc1","Lama4","Hspg2","Col4a2","Col1a2","Itga6","Itgb1","Itga1","App","Plxna2","Itgav","Itgb1","Sema6a","Ptprm","Pecam1","Jam2","Esam","Cdh5","APP","Ptprm","Cd74","Plxna4","Plxna2","Ptprm","Pecam1","Jam2","Esam","Cdh5","Cd74"

exprSet <- scRNAJe@assays$RNA@data
exprSet<-as.data.frame(t(exprSet))#转置

# 2. 计算某个基因和其它基因的相关性
y <- as.numeric(exprSet[,"Nucb2"])

col_pos <- match(Je_gene, colnames(exprSet))
col_pos=col_pos[!is.na(col_pos)]
df_existing <- exprSet[,col_pos]
colnames<-colnames(df_existing)
cor_data_df<- data.frame(colnames)
for(i in 1:length(colnames)){
  test_p<-cor.test(as.numeric(df_existing[,i]),y,method = "pearson")
  test_s<-cor.test(as.numeric(df_existing[,i]),y,method = "spearman")
  test_k<-cor.test(as.numeric(df_existing[,i]),y,method = "kendall")
  cor_data_df[i,2]<- test_s$estimate
  cor_data_df[i,3]<- test_s$p.value 
  cor_data_df[i,4]<- test_p$estimate
  cor_data_df[i,5]<- test_p$p.value
  cor_data_df[i,6]<- test_k$estimate
  cor_data_df[i,7]<- test_k$p.value 
}
names(cor_data_df)<-c("symbol","correlation-s","pvalue-s","correlation-p","pvalue-p","correlation-k","pvalue-k")
# cor_data_df %>% head()
write.csv(file='Je.Nucb2-interaction.cor3.csv',cor_data_df)


searchsabbr=c('SS','ECM-R','CCC')
# 比较分组的样本名称，用于保存文件
comparegroup=c('SeJe','SiJi','STJT')
# 比较分组的样本在object.list所在的位置进行分组，用于后续分析
groups=list(c(1,4),c(2,5),c(3,6))
searchs=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
for (i in 1:length(searchs)){
    search=searchs[i]
    abbr=searchsabbr[i]
    # search=searchs[1]
    # abbr=searchsabbr[1]
    # 整合Se Si ST Je Ji JT Secreted Signaling cellchat
    load(sprintf('Se_cellchat_%s.RData',search))
    Se_cellchat=cellchat
    load(sprintf('Si_cellchat_%s.RData',search))
    Si_cellchat=cellchat
    load(sprintf('ST_cellchat_%s.RData',search))
    ST_cellchat=cellchat
    load(sprintf('Je_cellchat_%s.RData',search))
    Je_cellchat=cellchat
    load(sprintf('Ji_cellchat_%s.RData',search))
    Ji_cellchat=cellchat
    load(sprintf('JT_cellchat_%s.RData',search))
    JT_cellchat=cellchat

    object.list <- list(Se=Se_cellchat,Je=Je_cellchat)
    #run netAnalysis_computeCentrality
    object.list<- lapply(object.list,function(x){
        x=netAnalysis_computeCentrality(x)})
    cellchat_m <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
    ####从宏观角度预测细胞通讯####
    #比较交互总数和交互强度
    compare=comparegroup[1]
    group=groups[i]
    # compare=comparegroup[3]
    # group=groups[[3]]
    pdf(file=sprintf("%sChatbubble_%s-AME.pdf",compare,abbr),height = 8, width = 8)
    par(mfrow = c(1,1), xpd=TRUE)
    netVisual_bubble(cellchat_m, sources.use =c('AP','myoFB','EC'),targets.use = c('AP','myoFB','EC'),comparison = c(1,2), angle.x = 45, font.size = 16) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 16))
    print(netVisual_bubble(cellchat_m, sources.use =c('AP','myoFB','EC'),targets.use = c('AP','myoFB','EC'),comparison =  c(1,2), angle.x = 45, font.size = 16)) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 16))
    dev.off()
}


load(sprintf('Se_cellchat_%s.RData','Secreted Signaling'))
Se_cellchat=cellchat
load(sprintf('Je_cellchat_%s.RData','Secreted Signaling'))
Je_cellchat=cellchat

object.list <- list(Se=Se_cellchat,Je=Je_cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
pathways.show <- c("VEGF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)

pdf(file=sprintf("SeJe-VEGF.pdf",compare,abbr),height = 8, width = 8) 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_aggregate(object.list[[1]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[1]))
netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[2], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[2]))
dev.off()

pathways.show <- c("IGF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)

pdf(file=sprintf("SeJe-IGF.pdf",compare,abbr),height = 8, width = 8) 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_aggregate(object.list[[1]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[1]))
netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[2], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[2]))
dev.off()


searchsabbr=c('SS','ECM-R','CCC')
comparegroup=c('SeJe','SiJi','STJT')
groups=list(c(1,4),c(2,5),c(3,6))
searchs=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
for (i in 1:length(searchs)){
    search=searchs[i]
    abbr=searchsabbr[i]
    load(sprintf('Se_cellchat_%s.RData',search))
    Se_cellchat=cellchat
    load(sprintf('Si_cellchat_%s.RData',search))
    Si_cellchat=cellchat
    load(sprintf('ST_cellchat_%s.RData',search))
    ST_cellchat=cellchat
    load(sprintf('Je_cellchat_%s.RData',search))
    Je_cellchat=cellchat
    load(sprintf('Ji_cellchat_%s.RData',search))
    Ji_cellchat=cellchat
    load(sprintf('JT_cellchat_%s.RData',search))
    JT_cellchat=cellchat

    object.list <- list(Je=Je_cellchat,Ji=Ji_cellchat,JT=JT_cellchat,Se=Se_cellchat,Si =Si_cellchat,ST=ST_cellchat)

    object.list<- lapply(object.list,function(x){
        x=netAnalysis_computeCentrality(x)})
    cellchat_m <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
    for (num in 1:length(groups)){
        compare=comparegroup[num]
        group=groups[[num]]
        pdf(file=sprintf("%sChatNumCell_%s.pdf",compare,abbr),height = 8, width = 8)
        par(mfrow = c(1,1), xpd=TRUE)
        netVisual_diffInteraction(cellchat_m, weight.scale = T,comparison = group)
        netVisual_diffInteraction(cellchat_m, weight.scale = T, measure = "weight",comparison = group)
        dev.off()
    }
}

