map <- read.table(
  file = "0828.group.txt",
  comment.char = "",
  row.names = 1,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t"
)

map2 <- map[id$sample,]
map2 <- as.data.frame(map2)
rownames(map2) <- id$sample
colnames(map2) <- c('Category1')
map2 <- na.omit(map2)
map2$Category1 <- factor(map2$Category1)

#####metabolite#####
metabolite <- read_excel("metabolite_data.xlsx")
ncol_metabolite <- ncol(metabolite)
metabolite <- metabolite[,c(2,17:ncol_metabolite)]

#####拆分pos/neg#####
mb_pos <- metabolite[c(metabolite$posneg=='pos'),]
mb_pos <- mb_pos[,c('name',id$id)]
rownames_mb_pos <- mb_pos$name
mb_pos <- mb_pos[,-1]
rownames(mb_pos) <- rownames_mb_pos
df_mb_pos <- t(mb_pos)


mb_neg <- metabolite[c(metabolite$posneg=='neg'),]
mb_neg <- mb_neg[,c('name',id$id)]
rownames_mb_neg <-mb_neg$name
mb_neg <- mb_neg[,-1]
rownames(mb_neg) <- rownames_mb_neg
df_mb_neg <- t(mb_neg)


####pos差异分析####
#if (!requireNamespace("BiocManager", quietly = TRUE))
# BiocManager::install("edgeR")

library(limma)
library(edgeR)

dge <- DGEList(counts=mb_pos)
group.list=map2$Category1
group.list=factor(group.list)
group.list
group.list=relevel(group.list,ref = "CONTROL")

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

DEG2=topTable(fit2, coef=1, adjust="BH")  
DEG = na.omit(DEG)   
write.csv(DEG,file = '0828.mb.P_C v.s. O.csv')

# pos_diff_mb <- DEG[DEG[, "P.Value"] < 0.05,]
# variables_pos <- rownames(pos_diff_mb)

# df_pos <- df_mb_pos[,c(variables_pos)]
# 筛选差异基因
## 我们使用|logFC| > 0.5，padj < 0.05（矫正后P值）
foldChange = 0.5
padj = 0.05
## 筛选出所有差异基因的结果
All_DEG <- DEG[(DEG$adj.P.Val < padj & (DEG$logFC>foldChange | DEG$logFC < (-foldChange))),]
#---------------------

## 我们发现竟然没有差异基因，这是应该我这边的数据是随机的结果，如果你的数据有这样的问题，你需要在仔细检查一下哦。
## 我们为了下面的操作正常进行，我们选用的P值（未矫正）进行筛选。
All_DEG <- DEG[(DEG$P.Value < padj & (DEG$logFC>foldChange | DEG$logFC < (-foldChange))),]
write.csv(All_DEG, "0828.mb.POS.all.DEG.csv")  ##输出差异基因数据集

# 筛选上调和下调的基因

diffup <-  All_DEG[(All_DEG$P.Value < padj & (All_DEG$logFC > foldChange)),]
write.csv(diffup, "0828.mb.POS.diffup.csv")
#
diffdown <- All_DEG[(All_DEG$P.Value < padj & (All_DEG < -foldChange)),]
write.csv(diffdown, "0828.mb.POS.diffdown.csv")
## 导入R包
library(ggplot2)
library(ggrepel)
##  绘制火山图
## 进行分类别
logFC <- DEG$logFC
deg.padj <- DEG$P.Value
data <- data.frame(logFC = logFC, padj = deg.padj)
data$group[(data$padj > padj | data$padj == "NA") | (data$logFC < foldChange) & data$logFC > -foldChange] <- "NS"
data$group[(data$padj <= padj & data$logFC > 0)] <-  "Up"
data$group[(data$padj <= padj & data$logFC < 0)] <- "Down"
x_lim <- max(logFC,-logFC)

# 开始绘图
pdf('0828.mb.POS.volcano.pdf',width = 7,height = 6.5)  ## 输出文件
label = subset(DEG,P.Value < padj & abs(logFC) > foldChange)
label1 = rownames(label)

colnames(DEG)[1] = 'log2FC'
Significant=ifelse((DEG$P.Value < padj & abs(DEG$log2FC)> foldChange), ifelse(DEG$log2FC > foldChange,"Up","Down"), "NS")

ggplot(DEG, aes(log2FC, -log10(P.Value)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#0072B5","grey","#BC3C28"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-0.5,0.5), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
 theme_bw()

dev.off()
## 导入R包
library(pheatmap)

## 提取差异基因的表达量
DEG_id < All_DEG  # 读取差异基因的文件
head(DEG_id)
## 匹配差异基因的表达量
DEG_id <- unique(DEG_id$X)
DEG_exp <- df[DEG_id,]
hmexp <- na.omit(DEG_exp)

## 样本注释信息 
annotation_col <- data.frame(Group = factor(map2$Category1))
rownames(annotation_col) <- colnames(hmexp)

##  绘制热图 
pdf(file = "0828.mb.POS.heatmap.pdf", height = 8, width = 12)
pheatmap(hmexp,
              annotation_col = annotation_col,
              color = colorRampPalette(c("green","black","red"))(50),
              cluster_cols = F,
              show_rownames = F,
              show_colnames = F,
              scale = "row", ## none, row, column
              fontsize = 12,
              fontsize_row = 12,
              fontsize_col = 6,
              border = FALSE)
print(p)
dev.off()

####neg差异分析#####

library(limma)
library(edgeR)

dge <- DGEList(counts=mb_neg)
group.list=map2$Category1
group.list=factor(group.list)
group.list=relevel(group.list,ref = "CONTROL")

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

DEG2=topTable(fit2, coef=1, adjust="BH")  
DEG = na.omit(DEG)   
write.csv(DEG,file = '0828.mb.N_C v.s. O.csv')

# neg_diff_mb <- DEG[DEG[, "P.Value"] < 0.05,]
# variables_neg <- rownames(neg_diff_mb)

# df_neg <- df_mb_neg[,c(variables_neg)]
# 筛选差异基因
## 我们使用|logFC| > 0.5，padj < 0.05（矫正后P值）

## 筛选出所有差异基因的结果
All_DEG <- DEG[(DEG$adj.P.Val < padj & (DEG$logFC>foldChange | DEG$logFC < (-foldChange))),]
#---------------------

## 我们发现竟然没有差异基因，这是应该我这边的数据是随机的结果，如果你的数据有这样的问题，你需要在仔细检查一下哦。
## 我们为了下面的操作正常进行，我们选用的P值（未矫正）进行筛选。
All_DEG <- DEG[(DEG$P.Value < padj & (DEG$logFC>foldChange | DEG$logFC < (-foldChange))),]
write.csv(All_DEG, "0828.mb.NEG.all.DEG.csv")  ##输出差异基因数据集

# 筛选上调和下调的基因

diffup <-  All_DEG[(All_DEG$P.Value < padj & (All_DEG$logFC > foldChange)),]
write.csv(diffup, "0828.mb.NEG.diffup.csv")
#
diffdown <- All_DEG[(All_DEG$P.Value < padj & (All_DEG < -foldChange)),]
write.csv(diffdown, "0828.mb.NEG.diffdown.csv")
## 导入R包
library(ggplot2)
library(ggrepel)
##  绘制火山图
## 进行分类别
logFC <- DEG$logFC
deg.padj <- DEG$P.Value
data <- data.frame(logFC = logFC, padj = deg.padj)
data$group[(data$padj > padj | data$padj == "NA") | (data$logFC < foldChange) & data$logFC > -foldChange] <- "NS"
data$group[(data$padj <= padj & data$logFC > 0)] <-  "Up"
data$group[(data$padj <= padj & data$logFC < 0)] <- "Down"
x_lim <- max(logFC,-logFC)

# 开始绘图
pdf('0828.mb.NEG.volcano.pdf',width = 7,height = 6.5)  ## 输出文件
label = subset(DEG,P.Value < padj & abs(logFC) > foldChange)
label1 = rownames(label)

colnames(DEG)[1] = 'log2FC'
Significant=ifelse((DEG$P.Value < padj & abs(DEG$log2FC)> foldChange), ifelse(DEG$log2FC > foldChange,"Up","Down"), "NS")

ggplot(DEG, aes(log2FC, -log10(P.Value)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#0072B5","grey","#BC3C28"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-0.5,0.5), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
 theme_bw()

dev.off()
## 导入R包
library(pheatmap)

## 提取差异基因的表达量
DEG_id < All_DEG  # 读取差异基因的文件
head(DEG_id)
## 匹配差异基因的表达量
DEG_id <- unique(DEG_id$X)
DEG_exp <- df[DEG_id,]
hmexp <- na.omit(DEG_exp)

## 样本注释信息 
annotation_col <- data.frame(Group = factor(map2$Category1))
rownames(annotation_col) <- colnames(hmexp)

##  绘制热图 
pdf(file = "0828.mb.NEG.heatmap.pdf", height = 8, width = 12)
pheatmap(hmexp,
              annotation_col = annotation_col,
              color = colorRampPalette(c("green","black","red"))(50),
              cluster_cols = F,
              show_rownames = F,
              show_colnames = F,
              scale = "row", ## none, row, column
              fontsize = 12,
              fontsize_row = 12,
              fontsize_col = 6,
              border = FALSE)
print(p)
dev.off()