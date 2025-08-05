rm(list=ls())


library(tidyr)
library(dplyr)
library(readxl) 
####联合分析的样本信息####
id <- read_excel("0828.id.xlsx")

####读取kraken2分析数据并处理#####
data <- read.csv(
  file = "./0828.S.count.tsv.all.tmp",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t",
  na.strings = "",
  fill=TRUE
)
####rawdata####
count_tmp <- data
taxonomy <- sapply(strsplit(count_tmp$taxonomy, "; g__"), function(x) x[length(x)])
taxonomy <- gsub("; s__", ";", taxonomy)
rownames(count_tmp) <- taxonomy
# 删除第一列OTU
count_tmp <- count_tmp[,-1]
# 删除最后一列taxonomy
df <- count_tmp %>% select(-taxonomy)

original_colnames <- colnames(df) 
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
colnames(df) <- cleaned_colnames  

#####剔除样本#####
# df <- df %>% select(-V350092088_96)
df_abun <- df
samples <- colnames(df)

df_abun <- df_abun[,c(id$sample)]


######提取Phylum#####
# otu <- data
# original_colnames <- colnames(otu) 
# cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# # 将修改后的列名赋值回数据框  
# colnames(otu) <- cleaned_colnames 
# otu <- otu %>% select(-taxonomy)
# rownames(otu) <- otu$`#OTU ID`
# otu <- otu[,-1]


# # 使用separate()函数分割combined_data列，并生成6列新数据  
# df_split <- separate(data, col = taxonomy, into = paste0("col", 1:7), sep = ";")  
# rownames(df_split) <- df_split$`#OTU ID`
# df_split <- df_split[,-1]
# last_7_cols <- df_split[, -(1:(ncol(df_split) - 7))] 
# colnames(last_7_cols) <- c('Kingdom','Phylum',
#                            'Class','Order',
#                            'Family','Genus',
#                            'Species')
# tax <- last_7_cols
# # if (!requireNamespace("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# # BiocManager::install("phyloseq")  
# 
# 
# 
# library(phyloseq)
# # otu_table <-read.csv("otu_table.csv",header = T,row.names = 1)
# # 
# # taxonomy_table <- read.csv("taxonomy_table.csv",header = T,row.names = 1)
# otu_table_phy <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
# taxonomy_table_phy <- tax_table(as.matrix(tax))    
# #合并
# physeq <- phyloseq(otu_table_phy, taxonomy_table_phy)   
# phylum_abundance <- tax_glom(physeq, taxrank = "Phylum")
# # 提取门的丰度表
# phylum_abundance_table <- otu_table(phylum_abundance)
# # 提取门的分类信息
# phylum_taxonomy <- tax_table(phylum_abundance)
# # 将丰度表和分类信息结合
# phylum_abundance_with_taxonomy <- cbind(as.data.frame(phylum_abundance_table), as.data.frame(phylum_taxonomy))
# # 保存结果到文件
# write.csv(phylum_abundance_with_taxonomy, file = "phylum_abundance_with_taxonomy.csv")
# 
# rownames(phylum_abundance_with_taxonomy) <- phylum_abundance_with_taxonomy$Phylum
# # 首先获取df的列数  
# ncol_df <- ncol(phylum_abundance_with_taxonomy)  
# # 然后选择除了最后7列之外的所有列  
# phylum_abundance_with_taxonomy <- phylum_abundance_with_taxonomy[, 1:(ncol_df - 7)]
# df_abun <- t(phylum_abundance_with_taxonomy)
# 
# # # # library(readxl) 
# id <- read_excel("id.xlsx")
# rownames(df_abun) <- id$id
# 
# df_abun <- df_abun[c(id2$id),]
# # group <- c(rep('C',4),rep('O',14))
# df1_row <- df_abun
# 
# df_abun <- t(df_abun)

####df_abun获取差异物种####
library("DESeq2")
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
# # 从样本信息中提取Category1信息
# category <- map["Category1"]
# # 移除NA值
# category <- na.omit(category)
# # category <- as.data.frame(group)
# rownames(category) <- 
# category <- category[c(id$id),]
# category <- as.data.frame(category)
# colnames(category) <- c('Category1')
# rownames(category) <- id2$id
# 在计数数据中加1，以避免在后续步骤中对零取对数
df_abun <- df_abun + 1
map2$Category1 <- factor(map2$Category1)

# 从计数数据和样本信息创建DESeqDataSet对象
dds <- DESeqDataSetFromMatrix(
  countData = df_abun,
  colData = map2,
  design = ~Category1 # 比较的分组方案
)
# 将负二项分布模型拟合到数据并进行差异表达测试
fit <- DESeq(dds, quiet = FALSE)
contrast = c("Category1", "OSAHS", "CONTROL")
# 提取差异表达分析的结果，将对比设置为将O组与C组进行比较
df <- results(
  fit,
  contrast = contrast
)

# 将结果转换为数据框架，并添加一个Features列，其中包含行名
df <- data.frame(
  Features = rownames(df),
  df,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
#####提取差异物种#####
pvalue <- 0.05
log2fc <- 2
pcol <- "padj"

df_up <- df[df$log2FoldChange > log2fc & df[, pcol] < pvalue,]
df_up <- df_up[complete.cases(df_up), ] 
df_upf <- df_up$Features
df_upf
df_down <- df[df$log2FoldChange < -log2fc & df[, pcol] < pvalue,]
df_down <- df_down[complete.cases(df_down), ] 
df_downf <- df_down$Features
df_downf
rbind_df <- rbind(df_up, df_down) 
variables <- rbind_df$Features


df1 <- t(df_abun)[,c(variables)]


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

pos_diff_mb <- DEG[DEG[, "P.Value"] < 0.05,]
variables_pos <- rownames(pos_diff_mb)

df_pos <- df_mb_pos[,c(variables_pos)]
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

neg_diff_mb <- DEG[DEG[, "P.Value"] < 0.05,]
variables_neg <- rownames(neg_diff_mb)

df_neg <- df_mb_neg[,c(variables_neg)]


######abun_pos相关分析#######
library(psych)  #psych包用于计算相关性、p值等信息
library(pheatmap) #用于绘制热图
library(reshape2) #reshape2包用于输出数据的整合处理

res <- corr.test(df1,df_pos,method = "spearman",alpha = 0.05) #method可选“pearson”、“spearman”、“kendall”
result_p <- res$p #提取p值
result_r <- res$r #提取cor值
p.out<-cbind(rownames(result_p),result_p) #输出p值
r.out<-cbind(rownames(result_r),result_r) #输出cor值
write.csv(p.out,"POS.pout.csv",row.names = F)
write.csv(r.out,"POS.rout.csv",row.names = F)

df_pr <-melt(result_r,value.name="cor")
df_pr$pvalue <-as.vector(result_p)  #宽格式转长格式
df_filt <- subset(df_pr,abs(df_pr$cor)>0.6&df_pr$pvalue<0.05)#筛选
# df_filt <- subset(df_pr,abs(df_pr$cor)>0.75)#筛选
# df_filt <- subset(df_pr,df_pr$pvalue<0.05)#筛选
write.csv(df_pr,"pos.melt_p_r.csv",row.names = F)#输出长格式
write.csv(df_filt,"pos.melt_p005_r06.csv",row.names = F)
# write.csv(df_filt,"melt_r075.csv",row.names = F)
# write.csv(df_filt,"melt_p005.csv",row.names = F)

if (!is.null(result_p)){
  ssmt <- result_p< 0.01
  result_p[ssmt] <-'**'
  smt <- result_p >0.01& result_p <0.05
  result_p[smt] <- '*'
  result_p[!ssmt&!smt]<- ''
} else {
  result_p <- F
}  #判断p值大小，若p<0.01，则'**'，若0.01<p<0.05,则'*'，否则无显著。
# pheatmap(result_r,file = "0828.abun_pos_cor.heatmap.pdf",width=16,height=8,scale = "row",cluster_row = T, cluster_col = T, border=NA,display_numbers = result_p, number_color = "black")#根据自己的喜好绘制热图

filt_r <- result_r[c(unique(df_filt$Var1)),c(unique(df_filt$Var2))]
# min(filt_r)

filt_p <- result_p[c(unique(df_filt$Var1)),c(unique(df_filt$Var2))]

pheatmap(filt_r,file = "0828.abun_pos_cor_filt.r6p005.heatmap.pdf",width=12,height=18,cluster_row = T, cluster_col = T, border=NA,display_numbers = filt_p, number_color = "black")#根据自己的喜好绘制热图

######abun_neg相关分析#######
library(psych)  #psych包用于计算相关性、p值等信息
library(pheatmap) #用于绘制热图
library(reshape2) #reshape2包用于输出数据的整合处理

res <- corr.test(df1,df_neg,method = "spearman",alpha = 0.05) #method可选“pearson”、“spearman”、“kendall”
result_p <- res$p #提取p值
result_r <- res$r #提取cor值
p.out<-cbind(rownames(result_p),result_p) #输出p值
r.out<-cbind(rownames(result_r),result_r) #输出cor值
write.csv(p.out,"neg.pout.csv",row.names = F)
write.csv(r.out,"neg.rout.csv",row.names = F)

df_pr <-melt(result_r,value.name="cor")
df_pr$pvalue <-as.vector(result_p)  #宽格式转长格式
df_filt <- subset(df_pr,abs(df_pr$cor)>0.6&df_pr$pvalue<0.05)#筛选
write.csv(df_pr,"neg.melt_p_r.csv",row.names = F)#输出长格式
write.csv(df_filt,"neg.melt_p005_r06.csv",row.names = F)

if (!is.null(result_p)){
  ssmt <- result_p< 0.01
  result_p[ssmt] <-'**'
  smt <- result_p >0.01& result_p <0.05
  result_p[smt] <- '*'
  result_p[!ssmt&!smt]<- ''
} else {
  result_p <- F
}  #判断p值大小，若p<0.01，则'**'，若0.01<p<0.05,则'*'，否则无显著。
filt_r <- result_r[c(unique(df_filt$Var1)),c(unique(df_filt$Var2))]
filt_p <- result_p[c(unique(df_filt$Var1)),c(unique(df_filt$Var2))]
pheatmap(filt_r,file = "0828.abun_neg_cor.r06p005.heatmap.pdf",width=12,height=18,cluster_row = T, cluster_col = T, border=NA,display_numbers = filt_p, number_color = "black")#根据自己的喜好绘制热图

####相关性网络pos####
##计算微生物类群丰度和代谢物鉴定强度的相关系数
library(Hmisc)
MAG<-df1
Enzyme<-df_pos

#计算群落组成与功能的相关性，以 spearman 相关系数为例
MAG_Enzyme_corr <- rcorr(as.matrix(MAG), as.matrix(Enzyme), type = 'spearman')

#相关系数 r 值和显著性 p 值矩阵
r <- MAG_Enzyme_corr$r
p <- MAG_Enzyme_corr$P

#只保留微生物丰度-代谢物的相关系数
#去除微生物-微生物、代谢物-代谢物之间的相关系数
r <- r[colnames(MAG),colnames(Enzyme)]
p <- p[colnames(MAG),colnames(Enzyme)]

#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
#该模式下，一定要注意负值的选择是否是合适的，因为很多情况下可能负相关无意义
r[abs(r) < 0.6] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
# p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
# p[p>=0.05] <- -1
# p[p<0.05 & p>=0] <- 1
# p[p==-1] <- 0
p[p>=0.05] <- 0
#根据上述筛选的 r 值和 p 值保留数据
z <- r * p

#再转换为对称矩阵，igraph 只能识别这种样式的邻接矩阵类型
z1 <- MAG_Enzyme_corr$r
z1[z1 != 0] <- 0
z1[rownames(z),colnames(z)] <- z
z1[colnames(z),rownames(z)] <- t(z)

#write.table(data.frame(z1, check.names = FALSE), 'MAG_Enzyme_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##获得网络
library(igraph)
#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物丰度和代谢物间的 spearman 相关系数
g <- graph.adjacency(z1, weighted = TRUE, mode = 'undirected')
g
#孤立节点的删除（删除度为 0 的节点）
g <- delete_vertices(g, names(degree(g)[degree(g) == 0]))
#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
#查看网络图
plot(g)
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write_graph(g, '0828.network_pos.graphml', format = 'graphml')


# dev.off()
####相关性网络neg####
##计算微生物类群丰度和代谢物鉴定强度的相关系数
library(Hmisc)
MAG<-df1
Enzyme<-df_neg

#计算群落组成与功能的相关性，以 spearman 相关系数为例
MAG_Enzyme_corr <- rcorr(as.matrix(MAG), as.matrix(Enzyme), type = 'spearman')

#相关系数 r 值和显著性 p 值矩阵
r <- MAG_Enzyme_corr$r
p <- MAG_Enzyme_corr$P

#只保留微生物丰度-代谢物的相关系数
#去除微生物-微生物、代谢物-代谢物之间的相关系数
r <- r[colnames(MAG),colnames(Enzyme)]
p <- p[colnames(MAG),colnames(Enzyme)]

#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
#该模式下，一定要注意负值的选择是否是合适的，因为很多情况下可能负相关无意义
r[abs(r) < 0.6] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
# p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
# p[p>=0.05] <- -1
# p[p<0.05 & p>=0] <- 1
# p[p==-1] <- 0
p[p>=0.05] <- 0
#根据上述筛选的 r 值和 p 值保留数据
z <- r * p

#再转换为对称矩阵，igraph 只能识别这种样式的邻接矩阵类型
z1 <- MAG_Enzyme_corr$r
z1[z1 != 0] <- 0
z1[rownames(z),colnames(z)] <- z
z1[colnames(z),rownames(z)] <- t(z)

#write.table(data.frame(z1, check.names = FALSE), 'MAG_Enzyme_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##获得网络
library(igraph)
#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物丰度和代谢物间的 spearman 相关系数
g <- graph.adjacency(z1, weighted = TRUE, mode = 'undirected')
g
#孤立节点的删除（删除度为 0 的节点）
g <- delete_vertices(g, names(degree(g)[degree(g) == 0]))
#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
#查看网络图

plot(g)
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write_graph(g, '0828.network_neg.graphml', format = 'graphml')

####pos_neg####
library(tibble)
df_pos <- as.data.frame(df_pos)
df_neg <- as.data.frame(df_neg)


# 转换行名为列，然后合并，最后（如果需要）将列转换回行名  
df_pos  <-df_pos  %>% rownames_to_column("RowName")  
df_neg <- df_neg %>% rownames_to_column("RowName")  
df_mb <- merge(df_pos, df_neg,by='RowName')  

# 如果需要将RowName转换回行名  
row.names(df_mb) <- df_mb$RowName  
df_mb$RowName <- NULL  # 删除RowName列
df_mb <- df_mb[,c('L-Theanine','Linatine')]

library(psych)  #psych包用于计算相关性、p值等信息
library(pheatmap) #用于绘制热图
library(reshape2) #reshape2包用于输出数据的整合处理

res <- corr.test(df1,df_mb,method = "spearman",alpha = 0.05) #method可选“pearson”、“spearman”、“kendall”
result_p <- res$p #提取p值
result_r <- res$r #提取cor值
p.out<-cbind(rownames(result_p),result_p) #输出p值
r.out<-cbind(rownames(result_r),result_r) #输出cor值
write.csv(p.out,"POS_NEG.all.pout.csv",row.names = F)
write.csv(r.out,"POS_NEG.all.rout.csv",row.names = F)

df_pr <-melt(result_r,value.name="cor")
df_pr$pvalue <-as.vector(result_p)  #宽格式转长格式
df_filt <- subset(df_pr,abs(df_pr$cor)>0.6&df_pr$pvalue<0.05)#筛选
# df_filt <- subset(df_pr,abs(df_pr$cor)>0.75)#筛选
# df_filt <- subset(df_pr,df_pr$pvalue<0.05)#筛选
write.csv(df_pr,"pos_neg.all.melt_p_r.csv",row.names = F)#输出长格式
write.csv(df_filt,"pos_neg.all.melt_p005_r06.csv",row.names = F)
# write.csv(df_filt,"melt_r075.csv",row.names = F)
# write.csv(df_filt,"melt_p005.csv",row.names = F)

if (!is.null(result_p)){
  ssmt <- result_p< 0.01
  result_p[ssmt] <-'**'
  smt <- result_p >0.01& result_p <0.05
  result_p[smt] <- '*'
  result_p[!ssmt&!smt]<- ''
} else {
  result_p <- F
}  #判断p值大小，若p<0.01，则'**'，若0.01<p<0.05,则'*'，否则无显著。
# pheatmap(result_r,file = "0828.abun_pos_cor.heatmap.pdf",width=16,height=8,scale = "row",cluster_row = T, cluster_col = T, border=NA,display_numbers = result_p, number_color = "black")#根据自己的喜好绘制热图

filt_r <- result_r[c(unique(df_filt$Var1)),c(unique(df_filt$Var2))]
# min(filt_r)

filt_p <- result_p[c(unique(df_filt$Var1)),c(unique(df_filt$Var2))]

pheatmap(filt_r,file = "0828.cor_ts.heatmap.pdf",width=6,height=8,cluster_row = T, cluster_col = T, border=NA,display_numbers = filt_p, number_color = "black")#根据自己的喜好绘制热图
####相关性网络pos_neg####
##计算微生物类群丰度和代谢物鉴定强度的相关系数
library(Hmisc)
MAG<-df1
Enzyme<-df_mb

#计算群落组成与功能的相关性，以 spearman 相关系数为例
MAG_Enzyme_corr <- rcorr(as.matrix(MAG), as.matrix(Enzyme), type = 'spearman')

#相关系数 r 值和显著性 p 值矩阵
r <- MAG_Enzyme_corr$r
p <- MAG_Enzyme_corr$P

#只保留微生物丰度-代谢物的相关系数
#去除微生物-微生物、代谢物-代谢物之间的相关系数
r <- r[colnames(MAG),colnames(Enzyme)]
p <- p[colnames(MAG),colnames(Enzyme)]

#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
#该模式下，一定要注意负值的选择是否是合适的，因为很多情况下可能负相关无意义
r[abs(r) < 0.6] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
# p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
# p[p>=0.05] <- -1
# p[p<0.05 & p>=0] <- 1
# p[p==-1] <- 0
p[p>=0.05] <- 0
#根据上述筛选的 r 值和 p 值保留数据
z <- r * p

#再转换为对称矩阵，igraph 只能识别这种样式的邻接矩阵类型
z1 <- MAG_Enzyme_corr$r
z1[z1 != 0] <- 0
z1[rownames(z),colnames(z)] <- z
z1[colnames(z),rownames(z)] <- t(z)

#write.table(data.frame(z1, check.names = FALSE), 'MAG_Enzyme_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##获得网络
library(igraph)
#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物丰度和代谢物间的 spearman 相关系数
g <- graph.adjacency(z1, weighted = TRUE, mode = 'undirected')
g
#孤立节点的删除（删除度为 0 的节点）
g <- delete_vertices(g, names(degree(g)[degree(g) == 0]))
#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
#查看网络图

plot(g)
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write_graph(g, '0828.network_cor_ts.graphml', format = 'graphml')

####O2PLS_pos####
# BiocManager::install("OmicsPLS")
library("OmicsPLS")
library(magrittr) # needs to be run every time you start R and want to use %>%
library(ggplot2)
df_abun_o <- df1_row
tax = scale(df_abun_o, scale=F)
mb_pos <- t(mb_pos)
met = scale(mb_pos, scale=F)

set.seed(123)
crossval_o2m(tax, met, 2:5,1:3,1:3,nr_folds = 10) #10折交叉验证
modelfit<-o2m(tax, met, 2, 3, 1)  #基于交叉验证结果确定成分数目参数
print (modelfit)

xj<- loadings(modelfit, "Xjoint", 1:2) %>% abs %>% rowSums
xj[-(order(xj,decreasing=T)[1:5])] = 0
xj <- sign(xj)
print(count(xj==1))
p1 <- plot(modelfit, loading_name="Xj", i=1, j=2, label = "c", use_ggplot2 = TRUE,
     alpha = xj,
     aes(label = stringr::str_sub(colnames(tax), start = 1)),size=4,col='red')+
  theme_bw() +
  coord_fixed(1, c(-1,1),c(-1,1)) +
  geom_point(alpha = 0.5+0.5*xj, col = 'blue',size=1.5) +
  labs(title = "taxonomy joint loadings",x = "First Joint Loadings", 
       y = "Second Joint Loadings") +
  theme(plot.title = element_text(face='bold'))
ggsave(
  plot = p1,
  "0828.O2PLS_pos_species.pdf",
  height = 8,
  width = 10
)        
yj<- loadings(modelfit, "Yjoint", 1:2) %>% abs %>% rowSums
yj[-(order(yj,decreasing=T)[1:5])] = 0
yj <- sign(yj)
print (yj)
p2 <- plot(modelfit, loading_name="Yj", i=1, j=2, label = "c", use_ggplot2 = TRUE,
     alpha = yj,
     aes(label = stringr::str_sub(colnames(met), start = 1)),size=4,col='red')+
  theme_bw() +
  coord_fixed(1, c(-1,1),c(-1,1)) +
  geom_point(alpha = 0.5+0.5*yj, col = 'blue',size=1.5) +
  labs(title ="metabolite joint loadings",
       x = "First Joint Loadings", y = "Second Joint Loadings") +
  theme(plot.title = element_text(face='bold'))      

# p1
# p2

ggsave(
  plot = p2,
  "0828.O2PLS_pos_metabolite.pdf",
  height = 8,
  width = 10
) 


####O2PLS_neg####
# BiocManager::install("OmicsPLS")
# library("OmicsPLS")
# library(magrittr) # needs to be run every time you start R and want to use %>%
# library(ggplot2)
# df_abun_o <- df1_row
# tax = scale(df_abun_o, scale=F)
mb_neg <- t(mb_neg)
met = scale(mb_neg, scale=F)

set.seed(123)
crossval_o2m(tax, met, 2:5,1:3,1:3,nr_folds = 10) #10折交叉验证
modelfit<-o2m(tax, met, 2, 3, 1)  #基于交叉验证结果确定成分数目参数
print (modelfit)

xj<- loadings(modelfit, "Xjoint", 1:2) %>% abs %>% rowSums
xj[-(order(xj,decreasing=T)[1:5])] = 0
xj <- sign(xj)
print(count(xj==1))
p3 <- plot(modelfit, loading_name="Xj", i=1, j=2, label = "c", use_ggplot2 = TRUE,
           alpha = xj,
           aes(label = stringr::str_sub(colnames(tax), start = 1)),size=4,col='red')+
  theme_bw() +
  coord_fixed(1, c(-1,1),c(-1,1)) +
  geom_point(alpha = 0.5+0.5*xj, col = 'blue',size=1.5) +
  labs(title = "taxonomy joint loadings",x = "First Joint Loadings", 
       y = "Second Joint Loadings") +
  theme(plot.title = element_text(face='bold'))
ggsave(
  plot = p3,
  "0828.O2PLS_neg_species.pdf",
  height = 8,
  width = 10
)        
yj<- loadings(modelfit, "Yjoint", 1:2) %>% abs %>% rowSums
yj[-(order(yj,decreasing=T)[1:5])] = 0
yj <- sign(yj)
print (yj)
p4 <- plot(modelfit, loading_name="Yj", i=1, j=2, label = "c", use_ggplot2 = TRUE,
           alpha = yj,
           aes(label = stringr::str_sub(colnames(met), start = 1)),size=4,col='red')+
  theme_bw() +
  coord_fixed(1, c(-1,1),c(-1,1)) +
  geom_point(alpha = 0.5+0.5*yj, col = 'blue',size=1.5) +
  labs(title ="metabolite joint loadings",
       x = "First Joint Loadings", y = "Second Joint Loadings") +
  theme(plot.title = element_text(face='bold'))      

# p1
# p2

ggsave(
  plot = p4,
  "0828.O2PLS_neg_metabolite.pdf",
  height = 8,
  width = 10
) 


####CCA####
library(vegan)
library(ggrepel)
library(ggplot2)

#冗余分析（redundancy analysis, RDA）或者典范对应分析（canonical correspondence analysis, CCA）,
#RDA将对应分析与多元回归分析相结合，每一步计算均与环境因子进行回归
#CCA是一种基于单峰模型的排序方法，而且在排序过程中结合多个环境因子
#RDA是基于线性模型，CCA是基于单峰模型，可以看哪种分析方法群落排序结果好进行选择
df_mb <- data.frame(df_mb_pos)
env <- data.frame(df1_row)

df_otu_cca <- cca(df_mb~., env)
cca <- df_otu_cca
ccascore <- scores(cca)

# #或者进行RDA分析
# RDA = rda(t(otu_tax[200:300,1:24]),design[,6:11], scale = TRUE)
# #获取RDA第一轴和第二轴,这里RDA不做可视化，以CCA为例子进行可视化
# RDA_sample<- scores(RDA,choices = 1:2, display = 'sp')
# RDA_env<-RDA$CCA$biplot[,1:2]

#获取CCA第一轴和第二轴
CCAE <- as.data.frame(cca$CCA$biplot[,1:2])
CCA1 <- ccascore$sites[,1]
CCA2 <- ccascore$sites[,2]
#构建样品第一轴和第二轴的数据框
plotdata <- data.frame(rownames(ccascore$sites), CCA1, CCA2)
colnames(plotdata) <- c("sample","CCA1","CCA2")
plotdata<-cbind(plotdata,env)
#获取第一轴和第二轴的解释度百分比，且保留一位小数
cca1 <- round(cca$CCA$eig[1]/sum(cca$CCA$eig)*100,1)
cca2 <- round(cca$CCA$eig[2]/sum(cca$CCA$eig)*100,1)
#绘制CAA分析结果
p_cca <- ggplot(plotdata,aes(x=CCA1,y=CCA2,color=plotdata$Treatment))+
  geom_point(size=4,aes())+
  stat_ellipse(aes(fill=plotdata$Treatment),geom = "polygon",size=0.6,level = 0.95,alpha = 0.1)+
  geom_segment(data=CCAE,aes(x = 0, y = 0, xend = CCAE[,1]*3.5, yend =  CCAE[,2]*3.5),
               arrow = arrow(length = unit(0.03, 'npc')),size =1,color="gray2")+
  geom_text(data = CCAE,aes(CCA1 * 4,CCA2 * 4,label = rownames(CCAE)),color = 'gray2',size = 4)+
  xlab(paste("CCA1 (",cca1,"%",")"))+ylab(paste("CCA2 (",cca2,"%",")"))+
  theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  scale_color_manual(values=c("#3FBDA7","#0172B6","#BD3C29","#F0965D"))+
  scale_fill_manual(values=c("#3FBDA7","#0172B6","#BD3C29","#F0965D"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_hline(aes(yintercept=0), colour="gray45",size=0.8, linetype="dashed")+
  geom_vline(aes(xintercept=0), colour="gray45",size=0.8, linetype="dashed")


p_cca


