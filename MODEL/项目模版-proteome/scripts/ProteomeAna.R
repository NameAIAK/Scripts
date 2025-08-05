rm(list=ls())
####load packages#####
library(readxl)
library(FactoMineR)
library(ggplot2)
library(ggrepel)
library(missMDA)
library(tidyr)
# install.packages('masstools')
# BiocManager::install("ComplexHeatmap")
library(scatterplot3d)
library(rgl)
library(EnhancedVolcano)
library(limma)
library(edgeR)
library(dplyr)
library(metpath)
# remotes:: install_gitlab("tidymass/metpath")
# remotes:: install_gitlab("tidymass/massdataset")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
options(timeout = 10000000000000000000000000)
# BiocManager::install("org.Mm.eg.db")
library(org.Hs.eg.db)
# library(org.Mm.eg.db) 
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggvenn)
library(VennDiagram)
library(stringr)

#先进行PCA聚类，选择样本后进行后续分析
norm_counts <- read.csv('data/lfq.proteins.csv')
colnames(norm_counts)[3] <- 'Accession_o'
norm_counts$Accession <- sapply(strsplit(norm_counts$Accession_o, "\\|"),function(x) x[1])
ACC_EN_GENE <- norm_counts[,c('Accession','Gene')]
norm_counts <- norm_counts[,c('Accession',
  'X1.Area',
  # 'X6.Area',
  'X13.Area',
  'X72.Area',
  'X75.Area',
  'X11.Area',
  'X17.Area',
  # 'X58.Area',
  'X70.Area',
  'X9.Area'
)]

####展示样本数据分组分布####
plot_data <- norm_counts
rownames(plot_data) <- plot_data$Accession
plot_data <- plot_data[,-1]
plot_data <- as.data.frame(t(plot_data))
plot_data$group <- factor(c(rep('F4',4),rep('F0',4)),levels = c('F4','F0'))
library(ggpubr)

# 设置每个组的颜色
group_colors <- c("F4" = "#BD3C29", "F0" = "#78D3AC")
# 首先设置比较的列表
compare_list <- list(
  c("F4","F0"))

# 使用ggplot2创建Shannon箱线图
p <- ggplot(plot_data, aes(x = group, y = Q16270, fill = group, colour = group)) +
  geom_boxplot(width = 0.5, alpha = 0.4, lwd = 1.15, outlier.shape = NA) +  # 调整箱的大小
  geom_jitter(width = 0.3, size = 3, alpha = 0.75) +  # 添加散点
  labs(y = "IGFBP7") +
  scale_fill_manual(values = group_colors) +  # 设置颜色
  scale_color_manual(values = group_colors) +
  theme_bw() +
  labs(x = NULL) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1.1),
        panel.grid.major = element_blank(), # 去除主网格线
        panel.grid.minor = element_blank()) +
  stat_compare_means(
    comparisons = compare_list,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE)# 添加检验结果
p
ggsave("IGFBP7.pdf", p, height = 3, width = 3)
# colnames(norm_counts) <- c('Accession','X1.Area','X6.Area','X13.Area','X72.Area',
#                            'X75.Area','X11.Area',
#                            'X17.Area','X58.Area','X70.Area','X9.Area')
accession <- norm_counts$Accession
norm_counts <- norm_counts[,-1]
rownames(norm_counts) <- accession
norm_counts[is.na(norm_counts)] <- 0

# PCAmat_t<-norm_counts

PCAmat_t<-t(norm_counts)
PCAmat_t[1:5,1:3]

#样本中基因表达值的 PCA 分析
gene.pca <- PCA(PCAmat_t, ncp = 3, scale.unit = FALSE, graph = FALSE)

#提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:3])
pca_sample$Sample=row.names(pca_sample)
#提取 PCA 前两轴的贡献度(22.46,19.13,12.62)
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2<- round(gene.pca$eig[2,2],2 )
pca_eig3<- round(gene.pca$eig[3,2],2 )
# head(pca_sample)
# rownames(pca_sample)
pca_sample$group<- factor(c(rep('F4',4),rep('F0',4)),levels = c('F4','F0'))
#pca_sample$group<- c(c(rep('B',5),rep('C',5),rep('M',5)))
str(pca_sample)
#####PCA12#####
p <- ggplot(data = pca_sample, aes(x =Dim.1 , y = Dim.2)) +
  geom_point(aes(color = group), size = 4.5) +  #根据样本坐标绘制二维散点图
   #自定义颜色
  scale_color_manual(values = c("#1E90FF","#E41A1C","green"))+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  #将 PCA 轴贡献度添加到坐标轴标题中

p
ggsave('PCA2Dpc12-nolabel.pdf',p)
library(ggrepel)
p <- p + 
  geom_text_repel(aes(label = Sample), size = 2.5, show.legend = FALSE, 
                  box.padding = unit(0.2, 'lines'))

p

# p <- p + 
#   stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)+ 
#   stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE) 

ggsave('PCA2Dpc12-label.pdf',p)

#####PCA13#####
p <- ggplot(data = pca_sample, aes(x =Dim.1 , y = Dim.3)) +
  geom_point(aes(color = group), size = 4.5) +  #根据样本坐标绘制二维散点图
  #自定义颜色
  scale_color_manual(values = c("#1E90FF","#E41A1C","green"))+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA3:', pca_eig3, '%'), color = '')  #将 PCA 轴贡献度添加到坐标轴标题中

p
ggsave('PCA2Dpc13-nolabel.pdf',p)
library(ggrepel)
p <- p + 
  geom_text_repel(aes(label = Sample), size = 2.5, show.legend = FALSE, 
                  box.padding = unit(0.2, 'lines'))

p

# p <- p + 
#   stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)+ 
#   stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE) 

ggsave('PCA2Dpc13-label.pdf',p)


#####PCA23#####
p <- ggplot(data = pca_sample, aes(x =Dim.2 , y = Dim.3)) +
  geom_point(aes(color = group), size = 4.5) +  #根据样本坐标绘制二维散点图
  #自定义颜色
  scale_color_manual(values = c("#1E90FF","#E41A1C","green"))+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PCA2:', pca_eig2, '%'), y = paste('PCA3:', pca_eig3, '%'), color = '')  #将 PCA 轴贡献度添加到坐标轴标题中

p
ggsave('PCA2Dpc23-nolabel.pdf',p)
library(ggrepel)
p <- p + 
  geom_text_repel(aes(label = Sample), size = 2.5, show.legend = FALSE, 
                  box.padding = unit(0.2, 'lines'))

p

# p <- p + 
#   stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)+ 
#   stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE) 

ggsave('PCA2Dpc23-label.pdf',p)

####3D
library(scatterplot3d)
dfGROUP<- data.frame(samples=colnames(norm_counts),group=colnames(norm_counts))
dfGROUP$group<-  factor(c(rep('F4',4),rep('F0',4)),levels = c('F4','F0'))
#dfGROUP$group<-  c(c(rep('B',5),rep('C',5),rep('M',5)))
# PCA计算
pca_result <- prcomp(t(norm_counts),
                     scale=F  # 一个逻辑值，指示在进行分析之前是否应该将变量缩放到具有单位方差
)
pca_result$x<-data.frame(pca_result$x)

colors <-factor(c(rep("#1E90FF",4),rep("#E41A1C",4)))
# colors <-c(rep("#1E90FF",3),rep("#E41A1C",3),rep('green',3))
#colors <-c(rep("#1E90FF",5),rep("#E41A1C",5),rep('green',5))


# 绘图
s3d <- scatterplot3d(pca_result$x[,1:3],
                     pch = 16,       # 点形状
                      color = colors,  # 点颜色
                     cex.symbols = 2 # 点大小
)

# 设置图例
legend("top",
       legend = c('F4','F0'),
       col =c("#1E90FF","#E41A1C") ,
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
PCAmat_t[1:3,1:3]

## pca分析
pca <- prcomp(t(PCAmat_t))

# 01.导入分类数据
dfGROUP<- data.frame(samples=rownames(PCAmat_t),group=rownames(PCAmat_t))
dfGROUP$group<- c(rep('F4',5),rep('F0',5))
#dfGROUP$group<- c(c(rep('B',5),rep('C',5),rep('M',5)))


## 
# 准备颜色
colors <-c(rep("#1E90FF",5),rep("#E41A1C",5))
#colors <-c(rep("#1E90FF",5),rep("#E41A1C",5),rep('green',5))
#ctrl,LPS,YJSTW16h,YJCTW16h,PJHASTW16h,PJHACTW16h,EZSTW16h,EZCTW16h

plot3d(pca_result$x[,1:3], # 取前三个主成分
       xlab="PC1", ylab="PC2", zlab="PC3", 
       col=colors, # 按groups填充颜色
       type="s", # 画球，'p' for points, 's' for spheres, 'l' for lines, 'h' for line segments 
       size=2, #球的大小
       lwd=2, box=T)
rgl.postscript("PCA3D_2.pdf",fmt="pdf",)


p <- ggplot(data = pca_sample, aes(x =Dim.1 , y = Dim.2)) +
  geom_point(aes(color = group), size = 4.5) +  #根据样本坐标绘制二维散点图
  #自定义颜色
  scale_color_manual(values = c("#1E90FF","#E41A1C","green"))+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA3:', pca_eig3, '%'), color = '')  #将 PCA 轴贡献度添加到坐标轴标题中

p

####差异分析####

P.Value <- 0.05
logFC <- 0


####M-C#####
norm_counts1 <- norm_counts
dge <- DGEList(counts=norm_counts1)
group.list=c(rep('F4',4),rep('F0',4))
group.list=factor(group.list)
group.list
group.list=relevel(group.list,ref = "F0")

design <- model.matrix(~0+group.list)
rownames(design)<-colnames(dge)
colnames(design)<-levels(group.list)
dge <- calcNormFactors(dge)   
v <- voom(dge,design, normalize="quantile")   
fit <- lmFit(v, design)    
constrasts = paste(rev(levels(group.list)),collapse = "-")
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 

fit2=contrasts.fit(fit,cont.matrix)   
fit2=eBayes(fit2)   
DEG = topTable(fit2, coef=constrasts,sort.by = "P", n=Inf)
# DEG_MC=topTable(fit2, coef=1, adjust="BH")  
DEG_MC = na.omit(DEG)   
write.csv(DEG_MC,file = 'proteome_F4vsF0.csv')
colnames(DEG_MC)
DEG_MC$change = ifelse(DEG_MC$P.Value< P.Value & abs(DEG_MC$logFC) > logFC, 
                   ifelse(DEG_MC$logFC> logFC ,'Up','Down'),
                   'Stable')
table(DEG_MC$change)
p <- ggplot(
  # 数据、映射、颜色
  DEG_MC, aes(x = logFC, y = -log10(P.Value), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="logFC",
       y="-log10(P.Value)")+
  scale_x_continuous(limits = c(-2.5, 2.5))+
  # scale_y_continuous(limits = c(0,30))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('F4 v.s. F0')
p
# 将火山图另存为SVG文件
ggsave(
  plot = p,
  "p_volcano.pdf",
  height = 4,
  width = 4
)

df_up <- DEG_MC[DEG_MC$logFC > logFC & DEG_MC[, 'P.Value'] < P.Value,]
df_up <- df_up[complete.cases(df_up), ] 
df_upf <- df_up$Accession
# df_upf

df_up <- df_up %>%  
  mutate(change = "UP")  
write.table(df_up, file = "proteins.up.csv", sep = ",", row.names = TRUE, quote = FALSE, na = "")

df_down <- DEG_MC[DEG_MC$logFC < -logFC & DEG_MC[, 'P.Value'] < P.Value,]
df_down <- df_down[complete.cases(df_down), ] 
df_downf <- df_down$Accession
# df_downf
df_down <- df_down %>%  
  mutate(change = "DOWN")
write.table(df_down, file = "proteins.down.csv", sep = ",", row.names = TRUE, quote = FALSE, na = "")

change <- bind_rows(df_up,df_down)
write.table(change, file = "proteins.change.csv", sep = ",", row.names = TRUE, quote = FALSE, na = "")

up_acc <- as.data.frame(df_upf)
colnames(up_acc) <- 'Accession'
up_acc <- merge(up_acc,ACC_EN_GENE,by='Accession')
up_acc <- up_acc[!(is.na(up_acc$`Gene`)),]
MC_UP <- up_acc
write.csv(MC_UP,file = 'Gene_UP.csv',row.names = FALSE)
# table(is.na(up_acc$`Ensembl Gene ID`))
# table(is.na(up_acc$`Gene Symbol`))
down_acc <- as.data.frame(df_downf)
colnames(down_acc) <- 'Accession'
down_acc <- merge(down_acc,ACC_EN_GENE,by='Accession')
down_acc <- down_acc[!(is.na(down_acc$`Gene`)),]
MC_DOWN <- down_acc
write.csv(MC_DOWN,file = 'Gene_DOWN.csv',row.names = FALSE)
####MC-enrichment####
####MC-UP####
gs = bitr(unique(up_acc$`Gene`), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")#org.Hs.eg.db;org.Mm.eg.db
####MC-GO_up####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_up <- data.frame(ego.bp@result)
ego.bp_up <- ego.bp_up[ego.bp_up$pvalue<0.05,]
# ego.bp_up <- ego.bp_up[c(1:10),]
# # ego.bp_up$GeneRatio_numeric <- sapply(ego.bp_up$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_up[order(ego.bp_up$GeneRatio_numeric),]$Description

ego.bp_up$Description<- factor(ego.bp_up$Description,levels =ego.bp_up[order(ego.bp_up$Count),]$Description)
write.csv(file="ego.bp_up_p0.05.csv",data.frame(ego.bp_up),row.names=F)

pdf(file="ego.bp_up_bar.pdf",width = 5,height = 5)
ggplot(ego.bp_up[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="GO Enrichment UP")+
  theme_bw()+ scale_y_discrete(labels=function(x) str_wrap(x, width=30))+theme(axis.text = element_text(size = 12))
dev.off()


####MC-kegg_up####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'hsa',pvalueCutoff = 0.05)#mmu;hsa
kk_up <- data.frame(kk@result)
kk_up <- kk_up[kk_up$pvalue<0.05,]
kk_up <- kk_up %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))
# kk_up <- ego.bp_up[c(1:10),]
# # ego.bp_up$GeneRatio_numeric <- sapply(ego.bp_up$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_up[order(ego.bp_up$GeneRatio_numeric),]$Description

kk_up$Description<- factor(kk_up$Description,levels =kk_up[order(kk_up$Count),]$Description)
write.csv(file="kk_up_p0.05.csv",data.frame(kk_up),row.names=F)

pdf(file="kk_up_bar.pdf",width = 5,height = 5)

ggplot(data = kk_up[1:9,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "KEGG_UP") +
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
gs = bitr(unique(down_acc$`Gene`), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
####MC-GO_down####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_down <- data.frame(ego.bp@result)
ego.bp_down <- ego.bp_down[ego.bp_down$pvalue<0.05,]
# ego.bp_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description
library(stringr)
ego.bp_down$Description<- factor(ego.bp_down$Description,levels =ego.bp_down[order(ego.bp_down$Count),]$Description)
write.csv(file="ego.bp_down_p0.05.csv",data.frame(ego.bp_down),row.names=F)

pdf(file="ego.bp_down_bar.pdf",width = 6,height =6)
ggplot(ego.bp_down[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="GO Enrichment DOWN")+
  theme_bw()+ scale_y_discrete(labels=function(x) str_wrap(x, width=30))+theme(axis.text = element_text(size = 12))
dev.off()


####MC-kegg_down####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'hsa',pvalueCutoff = 0.05)
kk_down <- data.frame(kk@result)
kk_down <- kk_down[kk_down$pvalue<0.05,]
kk_down <- kk_down %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))
# kk_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description

kk_down$Description<- factor(kk_down$Description,levels =kk_down[order(kk_down$Count),]$Description)
write.csv(file="kk_down_p0.05.csv",data.frame(kk_down),row.names=F)


pdf(file="kk_down_bar.pdf",width = 5,height = 5)

ggplot(data = kk_down[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "KEGG_down") +
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


####MC/TM inner gene####
MCUPinnerTMDOWN<- list('M vs C'=MC_UP$`Gene Symbol`,'T vs M'=TM_DOWN$`Gene Symbol`)
MCDOWNinnerTMUP<-list('M vs C'=MC_DOWN$`Gene Symbol`,'T vs M'=TM_UP$`Gene Symbol`)

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

pdf(file="MCUPinnerTMDOWN_GO.pdf",width = 7.5,height = 7)
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

pdf(file="MCDOWNinnerTMUP_GO.pdf",width =6.5,height = 7)
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


data1_deg <- read.csv('./proteins.change.csv')
colnames(data1_deg)[1] <- "Accession_N"
dataraw <- read.csv('data/lfq.proteins.csv')
dataraw19 <- dataraw[,1:19]
dataraw19$Accession_N <- sapply(strsplit(dataraw19$Accession, "\\|"),function(x) x[1])

df_M <- merge(dataraw19,data1_deg,by='Accession_N')
write.csv(file = 'ProteinsChange.csv',df_M)

########################
library(enrichplot)
library(clusterProfiler)
library(fgsea)
###获取ranks-geneList####
data <- merge(DEG_MC,ACC_EN_GENE,by='Accession')

cluster.genes<- data %>% arrange(desc(logFC)) %>% dplyr::select(Gene,logFC) #基因按logFC排序
ranks<- deframe(cluster.genes)
names(ranks)[names(ranks) == ""] <- "missing_gene"
ranks <- ranks[!duplicated(names(ranks))]
geneList = sort(ranks,decreasing = T)

###ref1##
genelist<- list.files('./gmt')
fgsea_sets<- list()
dir<- paste0('./gmt/',genelist[i])
temp<- read.gmt('./gmt/c5.all.v2024.1.Hs.symbols.gmt')

for (i in 1:length(names(table(temp$term))) ){
  # i=1
  termname <- names(table(temp$term))[i]
  fgsea_sets[[termname]]<- temp[temp$term==termname,]$gene
}
saveRDS(fgsea_sets,'C5.rda')
path <- as.data.frame(names(fgsea_sets))
colnames(path) <- 'name'
path$id <- rownames(path)
pathgrep <- path$name[grepl('FIBROSIS',path$name)]
# pathgrep <- path$name
pathgrepid <- path$id[grepl('FIBROSIS',path$name)]

geneset <- temp[temp$term %in% c(pathgrep,'GOBP_COMPLEMENT_ACTIVATION','HP_SYSTEMIC_LUPUS_ERYTHEMATOSUS','GOBP_BILE_ACID_SECRETION'), ]
# geneset <- temp[temp$term %in% c(pathgrep), ]
###ref2##
genelist<- list.files('./gmtfibrosis/')
fgsea_sets<- list()
for (i in 1:length(genelist)) {
  dir<- paste0('./gmtfibrosis/',genelist[i])
  temp<- read.gmt(dir)
  fgsea_sets[[names(table(temp$term))]]<- temp$gene
}
convert_to_geneset <- function(gene_sets) {
  # 对每个通路(term)展开其基因
  geneset_df <- data.frame(
    term = rep(names(gene_sets), lengths(gene_sets)),
    gene = unlist(gene_sets, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  
  # 转换为因子（如果需要）
  geneset_df$term <- as.factor(geneset_df$term)
  
  rownames(geneset_df) <- NULL
  return(geneset_df)
}
# 转换
geneSET <- convert_to_geneset(fgsea_sets)

####gsea富集####
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose = T,
             pAdjustMethod = 'BH',pvalueCutoff = 1,)

pdf(file = sprintf('gsea4rank%slabel.pdf',egmt@result$ID[1]),width = 6,height = 6)
gseaplot2(egmt,pvalue_table = T,
          geneSetID=1,#kegg共富集到54个通路，这里展示前3个通路
          subplots=1:3,
          ES_geom='line',title = egmt@result$ID[1]
          
)

dev.off()
pdf(file = 'gsea4dot.pdf',width = 6,height =6)
dotplot(egmt,split=".sign")+facet_grid(~.sign)
dev.off()
