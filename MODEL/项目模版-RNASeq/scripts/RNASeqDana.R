rm(list=ls())
####load packages#####
library(EnhancedVolcano)
library(limma)
library(edgeR)
library(readxl)
library(dplyr)
library(metPath)
library(ggplot2)
# library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggvenn)
library(VennDiagram)
library(org.Mm.eg.db) 
keytypes(org.Mm.eg.db)
library(ggpubr)
library(ggprism)
library(clusterProfiler)
library(pheatmap)
####合并所有样本数据####
dir=dir()
dat_all <- read.csv('./A1_2.dat.txt',sep = '\t')
# head(dat_all)
for (dat in dir[2:9]){
  print(dat)
  dat_content <- read.csv(dat,sep = '\t')
  dat_all <- merge(dat_all,dat_content,by='Geneid')
}
colnames(dat_all)
new_colnames <- gsub("X.data.snakemake.RNA.seq_up_mouse.bowtie2.", "", names(dat_all))
colnames(dat_all) <- new_colnames
new_colnames <- gsub(".bam", "", names(dat_all))
# 更新数据框的列名
colnames(dat_all) <- new_colnames

write.csv(dat_all,'1226_dat_all.csv',row.names = FALSE)

####读取所有样本数据####
dat_all <- read.csv('1226_dat_all.csv')
# 挑选样本
# dat_all <- dat_all[,c(1,2,3,5,7,8,10)]
gene<- dat_all$Geneid
for (i in 1:length(gene)) {
  temp<- strsplit(gene[i],'[.]')[[1]][1]
  gene[i]<- temp
}
data <- dat_all[,-1]
rownames(data)<- gene
gs<- bitr(rownames(data), fromType = "ENSEMBL", toType="SYMBOL",OrgDb = org.Mm.eg.db)
# table(duplicated(gs$ENSEMBL))
# table(duplicated(gs$SYMBOL))
####去除重复数据####
gs<- gs[!duplicated(gs$ENSEMBL),]
gs<- gs[!duplicated(gs$SYMBOL),]
data<- data[gs$ENSEMBL,]
rownames(data)<-gs$SYMBOL
###TMM标化
library(edgeR)
head(data)
for (i in 1:dim(data)[2]) {
  data[,i]<- as.numeric(data[,i])
}
group<- factor(c(rep('C',2),rep('M',2),rep('T',2)))
#group<- factor(c(rep('B',5),rep('C',5),rep('M',5)))
y<- DGEList(counts = data,group = group)
y<- calcNormFactors(y, method = "TMM")
y$samples
norm_counts<- cpm(y, normalized.lib.sizes = TRUE)
head(norm_counts)
write.csv(data,file = 'FTMMcounts.csv')
write.csv(norm_counts,file = 'TMM_norm_counts.csv')

####差异基因分析MC####
colnames(data)
#首先根据分组信息构建试验设计矩阵，分组信息中一定要是对照组在前，处理组在后
# data<- data[,-c(4,8,11)]
dd_CM<- data[,c(1:4)]
colnames(dd_CM)
group2<-  factor(rep(c('C', 'M'), each = 2))
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
write.table(lrt, 'MvsC_1.txt', sep = '\t', col.names = NA, quote = FALSE)

####差异基因分析TM####
colnames(data)
#首先根据分组信息构建试验设计矩阵，分组信息中一定要是对照组在前，处理组在后
# data<- data[,-c(4,8,11)]
dd_TM<- data[,c(3:6)]
colnames(dd_TM)
group2<-  factor(rep(c('M', 'T'), each = 2))
dgelist <- DGEList(counts =dd_TM, group = group2)
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
write.table(lrt, 'TvsM_1.txt', sep = '\t', col.names = NA, quote = FALSE)

# ####test####
# df2<- data.frame(1:9,1:9,1:9,1:9,1:9,1:9,1:9,1:9)
# colnames(df2)<- c('names','group','Tnf','Il6','Cxcl1','Il1b','Cxcl2','Nfe2l2')
# df2[,1]<- colnames(norm_counts)
# df2[,2]<- c(rep('C',3),rep('M',3),rep('T',3))
# df2[,3]<- as.numeric(norm_counts['Tnf',])
# df2[,4]<- as.numeric(norm_counts['Il6',])
# df2[,5]<- as.numeric(norm_counts['Cxcl1',])
# df2[,6]<- as.numeric(norm_counts['Il1b',])
# df2[,7]<- as.numeric(norm_counts['Cxcl2',])
# df2[,8]<- as.numeric(norm_counts['Nfe2l2',])
# 
# table(df2$group)
# df2$group<- factor(df2$group,levels = c('C','M','T'))
# p6 <- ggplot(df2,aes(x=group,y=Nfe2l2))+#指定数据
#   stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
#   geom_boxplot(aes(fill=group), #绘制箱线图函数
#                outlier.colour="white",size=0.8)+#异常点去除
#   theme(panel.background =element_blank(), #背景
#         axis.line=element_line(),#坐标轴的线设为显示
#         plot.title = element_text(size=14))+#图例位置
#   # scale_fill_manual(values=c("#ffc000","#a68dc8","blue"))+#指定颜色
#   geom_jitter(width = 0.2)+#添加抖动点
#   geom_signif(comparisons = list(c("M","C"),
#                                  c("T","M")),# 设置需要比较的组
#               map_signif_level = F, #是否使用星号显示
#               test = t.test, ##计算方法
#               
#               tip_length = c(c(0,0),
#                              c(0,0),
#                              c(0,0)),#横线下方的竖线设置
#               size=0.8,color="black")+
#   theme_prism(palette = "candy_bright",
#               base_fontface = "plain", # 字体样式，可选 bold, plain, italic
#               base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
#               base_size = 16,  # 图形的字体大小
#               base_line_size = 0.8, # 坐标轴的粗细
#               axis_text_angle = 45)+ # 可选值有 0，45，90，270
#   scale_fill_prism(palette = "candy_bright")+
#   theme(legend.position = 'none')#去除图例
# p6
# 
# ggarrange(p1,p2,p3,p4,p5,p6)
# 
# y_position = c(3,3.5,3.25)#图中横线位置 设置
# 
# 
# df2<- df2[-17,]
# df2<- df2[-c(1,7,15),]
# p6<- ggplot(df2,aes(x=group,y=Nfe2l2))+#指定数据
#   stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
#   geom_boxplot(aes(fill=group), #绘制箱线图函数
#                outlier.colour="white",size=0.8)+#异常点去除
#   theme(panel.background =element_blank(), #背景
#         axis.line=element_line(),#坐标轴的线设为显示
#         plot.title = element_text(size=14))+#图例位置
#   # scale_fill_manual(values=c("#ffc000","#a68dc8","blue"))+#指定颜色
#   geom_jitter(width = 0.2)+#添加抖动点
#   geom_signif(comparisons = list(c("A","M"),
#                                  c("B","M"),
#                                  c("C","M")),# 设置需要比较的组
#               map_signif_level = F, #是否使用星号显示
#               test = t.test, ##计算方法
#               y_position = c(28,27.5,27),
#               tip_length = c(c(0,0),
#                              c(0,0),
#                              c(0,0)),#横线下方的竖线设置
#               size=0.8,color="black")+
#   theme_prism(palette = "candy_bright",
#               base_fontface = "plain", # 字体样式，可选 bold, plain, italic
#               base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
#               base_size = 16,  # 图形的字体大小
#               base_line_size = 0.8, # 坐标轴的粗细
#               axis_text_angle = 45)+ # 可选值有 0，45，90，270
#   scale_fill_prism(palette = "candy_bright")+
#   theme(legend.position = 'none')#去除图例
# p6


####PCA####
PCAmat_t<-t(norm_counts)
PCAmat_t <- PCAmat_t[-c(3,5,8),]
# PCAmat_t <- PCAmat_t[ , which(apply(PCAmat_t, 2, var) != 0)]
# PCAmat_t[1:5,1:3]
library(FactoMineR)
library(ggplot2)
library(ggrepel)
#样本中基因表达值的 PCA 分析
gene.pca <- PCA(PCAmat_t, ncp = 3, scale.unit = FALSE, graph = FALSE)

#提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:3])
pca_sample$Sample=row.names(pca_sample)
#提取 PCA 前两轴的贡献度(22.46,19.13,12.62)
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )

head(pca_sample)
rownames(pca_sample)
pca_sample$group<- c(c(rep('C',2),rep('M',2),rep('T',2)))
#pca_sample$group<- c(c(rep('B',5),rep('C',5),rep('M',5)))

library(ggplot2)
p <- ggplot(data = pca_sample, aes(x =Dim.1 , y = Dim.2)) +
  geom_point(aes(color = group), size = 4.5) +  #根据样本坐标绘制二维散点图
  #自定义颜色
  scale_color_manual(values = c("#1E90FF","#E41A1C","green"))+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  #将 PCA 轴贡献度添加到坐标轴标题中

p

library(ggrepel)
p <- p + 
  geom_text_repel(aes(label = Sample), size = 2.5, show.legend = FALSE, 
                  box.padding = unit(0.2, 'lines'))

p

p <- p + stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)

p <- p + stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE) 
ggsave('PCA2D_2.pdf',p)
####3D
library(scatterplot3d)
# dfGROUP<- data.frame(samples=colnames(norm_counts),group=colnames(norm_counts))
# dfGROUP$group<-  c(c(rep('C',3),rep('M',3),rep('T',3)))
#dfGROUP$group<-  c(c(rep('B',5),rep('C',5),rep('M',5)))
# PCA计算
pca_result <- prcomp(PCAmat_t,
                     scale=F  # 一个逻辑值，指示在进行分析之前是否应该将变量缩放到具有单位方差
)
pca_result$x<-data.frame(pca_result$x)


colors <-c(rep("#1E90FF",2),rep("#E41A1C",2),rep('green',2))
#colors <-c(rep("#1E90FF",5),rep("#E41A1C",5),rep('green',5))


# 绘图
s3d <- scatterplot3d(pca_result$x[,1:3],
                     pch = 16,       # 点形状
                     color = colors,  # 点颜色
                     cex.symbols = 2 # 点大小
)
# 设置图例
legend("top",
       legend = c('C','M','T'),
       col =c("#1E90FF","#E41A1C","green") ,
       pch = 16,
       inset = -0.1,
       xpd = TRUE,
       horiz = TRUE)
#legend("top",legend = c('B','C','M'),col =c("#1E90FF","#E41A1C","green") ,pch = 16,inset = -0.1,xpd = TRUE,horiz = TRUE)
# 设置文字标注
text(s3d$xyz.convert(pca_result$x[,c(1,2,3)]),
     labels = row.names(pca_result$x),
     cex = 0.6,col = "black")


####3D_2
# 导入包
library(rgl)
options(stringsAsFactors = FALSE) #禁止chr转成factor
# PCAmat_t[1:3,1:3]

## pca分析
# pca <- prcomp(t(PCAmat_t))

# 01.导入分类数据
# dfGROUP<- data.frame(samples=rownames(PCAmat_t),group=rownames(PCAmat_t))
# dfGROUP$group<- c(c(rep('C',3),rep('M',3),rep('T',3)))
#dfGROUP$group<- c(c(rep('B',5),rep('C',5),rep('M',5)))


## 
# 准备颜色
colors <-c(rep("#1E90FF",3),rep("#E41A1C",3),rep('green',3))
#colors <-c(rep("#1E90FF",5),rep("#E41A1C",5),rep('green',5))
#ctrl,LPS,YJSTW16h,YJCTW16h,PJHASTW16h,PJHACTW16h,EZSTW16h,EZCTW16h

plot3d(pca_result$x[,1:3], # 取前三个主成分
       xlab="PC1", ylab="PC2", zlab="PC3", 
       col=colors, # 按groups填充颜色
       type="s", # 画球，'p' for points, 's' for spheres, 'l' for lines, 'h' for line segments 
       size=2, #球的大小
       lwd=2, box=T)

#####################3d_3
library(pca3d)
data(metabo)

rotation<- pca_result$rotation[,1:3]
colnames(rotation)
PC1<- rotation[,'PC1']
names(PC1)<- rownames(rotation)
PC1[1:5]
sort(PC1,decreasing = T)[1:10]

##
PC2<- rotation[,'PC2']
names(PC2)<- rownames(rotation)
PC2[1:5]
sort(PC2,decreasing = T)[1:10]

##
PC3<- rotation[,'PC3']
names(PC3)<- rownames(rotation)
PC3[1:5]
sort(PC3,decreasing = T)[1:10]


####火山图MC####
library(ggplot2)

MC<- read.table('MvsC_1.txt')
colnames(MC)
MC$gene<- rownames(MC)

# 设置p_value和logFC的阈值
cut_off_FDR = 0.05  #统计显著性
cut_off_PValue= 0.05
cut_off_logFC = 0         #差异倍数值
MC$change<- 'Stable'
# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
MC$change = ifelse(MC$PValue< cut_off_PValue & abs(MC$logFC) > cut_off_logFC, 
                   ifelse(MC$logFC> cut_off_logFC ,'Up','Down'),
                   'Stable')
head(MC)
table(MC$change)
write.csv(file = 'MC_DEG.csv',MC)
p <- ggplot(
  # 数据、映射、颜色
  MC, aes(x = logFC, y = -log10(PValue), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_PValue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (PValue)")+
  scale_x_continuous(limits = c(-5, 5))+
  scale_y_continuous(limits = c(0,10))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('M v.s. C')
p

# 将需要标记的基因放置在label列(logFC >= 5)
# library(ggrepel)
# MC$label <- ifelse(MC$PValue < cut_off_FDR & abs(MC$logFC) >= cut_off_logFC,
#                    as.character(MC$gene), "")
# 
# 
# p <- p + geom_label_repel(data = MC, aes(x = MC$logFC, 
#                                     y = -log10(MC$PValue), 
#                                     label = label),
#                      size = 3, box.padding = unit(0.5, "lines"),
#                      point.padding = unit(0.8, "lines"), 
#                      segment.color = "black", 
#                      show.legend = FALSE)
ggsave(
  plot = p,
  "MC.p_volcano.pdf",
  height = 6,
  width = 6
)
####火山图TM####
library(ggplot2)

TM<- read.table('TvsM_1.txt')
colnames(TM)
TM$gene<- rownames(TM)

TM$change<- 'Stable'
# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
TM$change = ifelse(TM$PValue< cut_off_PValue & abs(TM$logFC) > cut_off_logFC, 
                   ifelse(TM$logFC> cut_off_logFC ,'Up','Down'),
                   'Stable')
head(TM)
table(TM$change)
write.csv(file = 'TM_DEG.csv',TM)
p <- ggplot(
  # 数据、映射、颜色
  TM, aes(x = logFC, y = -log10(PValue), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_PValue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (PValue)")+
  scale_x_continuous(limits = c(-5, 5))+
  scale_y_continuous(limits = c(0,10))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('T v.s. M')
p

# 将需要标记的基因放置在label列(logFC >= 5)
# library(ggrepel)
# TM$label <- ifelse(TM$PValue < cut_off_FDR & abs(TM$logFC) >= cut_off_logFC,
#                    as.character(TM$gene), "")
# 
# 
# p <- p + geom_label_repel(data = TM, aes(x = TM$logFC, 
#                                     y = -log10(TM$PValue), 
#                                     label = label),
#                      size = 3, box.padding = unit(0.5, "lines"),
#                      point.padding = unit(0.8, "lines"), 
#                      segment.color = "black", 
#                      show.legend = FALSE)
ggsave(
  plot = p,
  "TM.p_volcano.pdf",
  height = 6,
  width = 6
)

####热图MC####
# dd<- data[c(MC[MC$change=='Up',]$gene,MC[MC$change=='Down',]$gene),1:6]
# 
# pheatmap(
#   dd,  # 丰度数据
#   scale='column',
#   file = "MC.heatmap.pdf",  # 输出文件名
#   # file = "0828.all.heatmap.pdf",
#   fontsize = 10,  # 字体大小
#   border_color = "white",  # 边框颜色
#   # color = colorRampPalette(  # 颜色渐变
#   # colors = c("#20A4B2", "#F9F7F7", "#FF0000")
#   # )(100),
#   # cluster_cols =TRUE,  # 是否对列进行聚类
#   cluster_cols =FALSE,
#   clustering_distance_cols = "euclidean",  # 列聚 类距离计算方法
#   # cluster_rows = TRUE,  # 是否对行进行聚类
#   cluster_rows = FALSE,
#   # scale = 'column',
#   cellwidth = 10,  # 单元格宽度
#   cellheight = 10,  # 单元格高度
#   # width = 8, height = 8,
#   show_rownames = TRUE,  # 是否显示行名
#   show_colnames = TRUE,  # 是否显示列名
#   # annotation_col  = group,  # 行注释
#   gaps_row = c(length(MC[MC$change=='Up',]$gene))
# )
# dev.off()
####热图TM####
####MC-enrichment####
####MC-UP####
MC_UP <- MC[MC$change=='Up',]
gs = bitr(unique(MC_UP$gene), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
####MC-GO_up####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_up <- data.frame(ego.bp@result)
ego.bp_up <- ego.bp_up[ego.bp_up$pvalue<0.05,]
# ego.bp_up <- ego.bp_up[c(1:10),]
# # ego.bp_up$GeneRatio_numeric <- sapply(ego.bp_up$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_up[order(ego.bp_up$GeneRatio_numeric),]$Description

ego.bp_up$Description<- factor(ego.bp_up$Description,levels =ego.bp_up[order(ego.bp_up$Count),]$Description)
write.csv(file="MC_ego.bp_up_p0.05.csv",data.frame(ego.bp_up),row.names=F)

pdf(file="MC_ego.bp_up_bar.pdf",width = 6,height = 7)
ggplot(ego.bp_up[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="MC GO Enrichment UP")+
  theme_bw()
dev.off()


####MC-kegg_up####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'mmu',pvalueCutoff = 0.05)
kk_up <- data.frame(kk@result)
kk_up <- kk_up[kk_up$pvalue<0.05,]
kk_up <- kk_up %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))
# kk_up <- ego.bp_up[c(1:10),]
# # ego.bp_up$GeneRatio_numeric <- sapply(ego.bp_up$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_up[order(ego.bp_up$GeneRatio_numeric),]$Description

kk_up$Description<- factor(kk_up$Description,levels =kk_up[order(kk_up$Count),]$Description)
write.csv(file="MC_kk_up_p0.05.csv",data.frame(kk_up),row.names=F)

pdf(file="MC_kk_up_bar.pdf",width = 7,height = 7)

ggplot(data = kk_up[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "MC_KEGG_UP") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()

####MC-DOWN####
MC_DOWN <- MC[MC$change=='Down',]
gs = bitr(unique(MC_DOWN$gene), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
####MC-GO_down####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_down <- data.frame(ego.bp@result)
ego.bp_down <- ego.bp_down[ego.bp_down$pvalue<0.05,]
# ego.bp_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description

ego.bp_down$Description<- factor(ego.bp_down$Description,levels =ego.bp_down[order(ego.bp_down$Count),]$Description)
write.csv(file="MC_ego.bp_down_p0.05.csv",data.frame(ego.bp_down),row.names=F)

pdf(file="MC_ego.bp_down_bar.pdf",width = 8.5,height = 7)
ggplot(ego.bp_down[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="MC GO Enrichment DOWN")+
  theme_bw()
dev.off()


####MC-kegg_down####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'mmu',pvalueCutoff = 0.05)
kk_down <- data.frame(kk@result)
kk_down <- kk_down[kk_down$pvalue<0.05,]
kk_down <- kk_down %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))
# kk_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description

kk_down$Description<- factor(kk_down$Description,levels =kk_down[order(kk_down$Count),]$Description)
write.csv(file="MC_kk_down_p0.05.csv",data.frame(kk_down),row.names=F)


pdf(file="MC_kk_down_bar.pdf",width = 7,height = 7)

ggplot(data = kk_down[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "MC_KEGG_down") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()
####TM-enrichment####
####TM-UP####
TM_UP <- TM[TM$change=='Up',]
gs = bitr(unique(TM_UP$gene), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
####TM-GO_up####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_up <- data.frame(ego.bp@result)
ego.bp_up <- ego.bp_up[ego.bp_up$pvalue<0.05,]
# ego.bp_up <- ego.bp_up[c(1:10),]
# # ego.bp_up$GeneRatio_numeric <- sapply(ego.bp_up$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_up[order(ego.bp_up$GeneRatio_numeric),]$Description

ego.bp_up$Description<- factor(ego.bp_up$Description,levels =ego.bp_up[order(ego.bp_up$Count),]$Description)
write.csv(file="TM_ego.bp_up_p0.05.csv",data.frame(ego.bp_up),row.names=F)

pdf(file="TM_ego.bp_up_bar.pdf",width = 6,height = 7)
ggplot(ego.bp_up[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="TM GO Enrichment UP")+
  theme_bw()
dev.off()


####TM-kegg_up####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'mmu',pvalueCutoff = 0.05)
kk_up <- data.frame(kk@result)
kk_up <- kk_up[kk_up$pvalue<0.05,]
kk_up <- kk_up %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))

# kk_up <- ego.bp_up[c(1:10),]
# # ego.bp_up$GeneRatio_numeric <- sapply(ego.bp_up$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_up[order(ego.bp_up$GeneRatio_numeric),]$Description

kk_up$Description<- factor(kk_up$Description,levels =kk_up[order(kk_up$Count),]$Description)
write.csv(file="TM_kk_up_p0.05.csv",data.frame(kk_up),row.names=F)


pdf(file="TM_kk_up_bar.pdf",width = 7,height = 7)

ggplot(data = kk_up[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "TM_KEGG_UP") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()

####TM-DOWN####
TM_DOWN <- TM[TM$change=='Down',]
gs = bitr(unique(TM_DOWN$gene), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
####TM-GO_down#####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_down <- data.frame(ego.bp@result)
ego.bp_down <- ego.bp_down[ego.bp_down$pvalue<0.05,]
# ego.bp_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description

ego.bp_down$Description<- factor(ego.bp_down$Description,levels =ego.bp_down[order(ego.bp_down$Count),]$Description)
write.csv(file="TM_ego.bp_down_p0.05.csv",data.frame(ego.bp_down),row.names=F)

pdf(file="TM_ego.bp_down_bar.pdf",width = 6,height = 7)
ggplot(ego.bp_down[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="TM GO Enrichment DOWN")+
  theme_bw()
dev.off()


####TM-kegg_down####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'mmu',pvalueCutoff = 0.05)
kk_down <- data.frame(kk@result)
kk_down <- kk_down[kk_down$pvalue<0.05,]
kk_down <- kk_down %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))
# kk_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description

kk_down$Description<- factor(kk_down$Description,levels =kk_down[order(kk_down$Count),]$Description)
write.csv(file="TM_kk_down_p0.05.csv",data.frame(kk_down),row.names=F)


pdf(file="TM_kk_down_bar.pdf",width = 7,height = 7)

ggplot(data = kk_down[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "TM_KEGG_down") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()


####venn####
library(ggvenn)
####MC/TM inner gene####
MCUPinnerTMDOWN<- list('M vs C'=MC_UP$gene,'T vs M'=TM_DOWN$gene)
MCDOWNinnerTMUP<-list('M vs C'=MC_DOWN$gene,'T vs M'=TM_UP$gene)

####MCUPinnerTMDOWN####
p_MCUPinnerTMDOWN <- ggvenn(MCUPinnerTMDOWN,c('M vs C', 'T vs M'),show_percentage = T,
                            stroke_color = "white",
                            fill_color = c("#BC3C29FF","#0072B5FF","#E18727FF","#FFDC91FF",
                                           "#EE4C97FF","#20854EFF","#7876B1FF","#6F99ADFF" ),
                            set_name_color =c("#E41A1C","#1E90FF"))
ggsave(p_MCUPinnerTMDOWN,filename = 'MCUPinnerTMDOWN.venn.pdf',width = 6,height = 4)

inner1 <- get.venn.partitions(MCUPinnerTMDOWN)
inner1_gene <- data.frame(inner1[[1,'..values..']])
colnames(inner1_gene) <- 'MCUPinnerTMDOWN_gene'
write.table(inner1_gene, 'MCUPinnerTMDOWN_gene_inter.csv', row.names = FALSE, col.names = TRUE,sep = ',', quote = FALSE)

gs = bitr(inner1_gene$MCUPinnerTMDOWN_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
####inner1GO#####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_inner1 <- data.frame(ego.bp@result)
ego.bp_inner1 <- ego.bp_inner1[ego.bp_inner1$pvalue<0.05,]
# ego.bp_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description
ego.bp_inner1$Description<- factor(ego.bp_inner1$Description,levels =ego.bp_inner1[order(ego.bp_inner1$Count),]$Description)
write.csv(file="MCUPinnerTMDOWN_GO_p0.05.csv",data.frame(ego.bp_inner1),row.names=F)

pdf(file="MCUPinnerTMDOWN_GO.pdf",width = 6,height = 7)
ggplot(ego.bp_inner1[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="MCUPinnerTMDOWN_GO")+
  theme_bw()
dev.off()

####inner1KEGG####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'mmu',pvalueCutoff = 0.05)
kk_inner1 <- data.frame(kk@result)
kk_inner1 <- kk_inner1[kk_inner1$pvalue<0.05,]
kk_inner1 <- kk_inner1 %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))
# kk_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description

kk_inner1$Description<- factor(kk_inner1$Description,levels =kk_inner1[order(kk_inner1$Count),]$Description)
write.csv(file="MCUPinnerTMDOWN_KEGG_p0.05.csv",data.frame(kk_inner1),row.names=F)


pdf(file="MCUPinnerTMDOWN_KEGG_bar.pdf",width = 7,height = 7)

ggplot(data = kk_inner1[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "MCUPinnerTMDOWN_KEGG") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()

####MCDOWNinnerTMUP####
p_MCDOWNinnerTMUP <- ggvenn(MCDOWNinnerTMUP,c('M vs C', 'T vs M'),show_percentage = T,
                            stroke_color = "white",
                            fill_color = c("#BC3C29FF","#0072B5FF","#E18727FF","#FFDC91FF",
                                           "#EE4C97FF","#20854EFF","#7876B1FF","#6F99ADFF" ),
                            set_name_color =c("#E41A1C","#1E90FF"))
ggsave(p_MCDOWNinnerTMUP,filename = 'MCDOWNinnerTMUP.venn.pdf',width = 6,height = 4)

inner2 <- get.venn.partitions(MCDOWNinnerTMUP)
inner2_gene <- data.frame(inner2[[1,'..values..']])
colnames(inner2_gene) <- 'MCDOWNinnerTMUP_gene'
write.table(inner2_gene, 'MCDOWNinnerTMUP_gene_inter.csv', row.names = FALSE, col.names = TRUE,sep = ',', quote = FALSE)

gs = bitr(inner2_gene$MCDOWNinnerTMUP_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
####inner2GO#####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_inner2 <- data.frame(ego.bp@result)
ego.bp_inner2 <- ego.bp_inner2[ego.bp_inner2$pvalue<0.05,]
# ego.bp_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description
ego.bp_inner2$Description<- factor(ego.bp_inner2$Description,levels =ego.bp_inner2[order(ego.bp_inner2$Count),]$Description)
write.csv(file="MCDOWNinnerTMUP_GO_p0.05.csv",data.frame(ego.bp_inner2),row.names=F)

pdf(file="MCDOWNinnerTMUP_GO.pdf",width = 6,height = 7)
ggplot(ego.bp_inner2[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="MCDOWNinnerTMUP_GO")+
  theme_bw()
dev.off()

####inner2KEGG####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'mmu',pvalueCutoff = 0.05)
kk_inner2 <- data.frame(kk@result)
kk_inner2 <- kk_inner2[kk_inner2$pvalue<0.05,]
kk_inner2 <- kk_inner2 %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))
# kk_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description

kk_inner2$Description<- factor(kk_inner2$Description,levels =kk_inner2[order(kk_inner2$Count),]$Description)
write.csv(file="MCDOWNinnerTMUP_KEGG_p0.05.csv",data.frame(kk_inner2),row.names=F)

pdf(file="MCDOWNinnerTMUP_KEGG_bar.pdf",width = 7,height = 7)
ggplot(data = kk_inner2[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "MCDOWNinnerTMUP_KEGG") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()







