#!/usr/bin/Rscript
####load packages#####
library(limma)
library(edgeR)
library(readxl)
library(dplyr)
library(ggplot2)
library(DOSE)
library(enrichplot)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(ggpubr)
library(ggprism)
library(clusterProfiler)
library(writexl)

# cmd:Rscript ../RNAseqDEGS.R GSE188774 Hs Control Treat 3 2 3 5 ../01_Rawdata/GSE188774_readcount_genename.xls Pvalue 0 0.05
# args[1] 项目名称：GSE188774
# args[2] 数据库名称缩写：Hs/Mm
# args[3] 对照组名称：Control
# args[4] 治疗组/实验组名称：Treat
# args[5] Control样本数：3
# args[6] Control开始数据列：2
# args[7] Treat样本数：3
# args[8] Treat开始数据列：5
# args[9] count文件：../01_Rawdata/GSE188774_readcount_genename.xls
# args[10]  显著性采用的值Pvalue/FDR
# args[11]  FC筛选值：0
# args[12]  P筛选值：0.05

# 查看操作路径
getwd()
# 获取参数
args <- commandArgs(trailingOnly = TRUE)
projectname <- args[1]#GSE188774
# 数据库为Hs/Mm
db_short <- args[2]#Hs
if (db_short=='Hs'){
  db='org.Hs.eg.db'
}else if (db_short=='Mm'){
  db='org.Mm.eg.db'
}
# 对照组名称
Controlname <-args[3] #'Control'
# 治疗组或实验组名称
Treatname <- args[4] #'Treat'
# 对照组样本数量
Controlsamplenum <- as.numeric(args[5]) #3
Controlsamplestartcol <- as.numeric(args[6])#2
Controlsampleendcol <- as.numeric(Controlsamplestartcol)+as.numeric(Controlsamplenum)-1
# 治疗组样本数量
Treatsamplenum <- as.numeric(args[7]) #3
Treatsamplestartcol <- as.numeric(args[8]) #4
Treatsampleendcol <- as.numeric(Treatsamplestartcol)+as.numeric(Treatsamplenum)-1
# count路径-文件名
countfile <- args[9] #'count.xls'
# 获取count文件类型
countfiletype <- tools::file_ext(countfile)
# 读取数据
if (countfiletype=='xls'){
  count <- read.csv2(countfile,sep = '\t')
} else if (file_extension %in% c("xlsx")) {
  print("这是Excel文件")
} else {
  print("未知文件类型")
}
P <- args[10]# Pvalue
# P <- 'FDR'
cut_off_logFC = as.numeric(args[11]) #0         #差异倍数值

# 挑选样本
dat_all <- count[,c(1,Controlsamplestartcol:Controlsampleendcol,Treatsamplestartcol:Treatsampleendcol)]
gene<- dat_all$gene_id
for (i in 1:length(gene)) {
  temp<- strsplit(gene[i],'[.]')[[1]][1]
  gene[i]<- temp
}
write.csv(dat_all,file = sprintf('%s_dat_all_counts.csv',projectname))
# 删除gene_id列
data <- dat_all[,-1]
rownames(data)<- gene
# gene_id to gene_name
gs<- bitr(rownames(data), fromType = "ENSEMBL", toType="SYMBOL",OrgDb = db)
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
# 创建分组
group<- factor(c(rep(Controlname,Controlsamplenum),rep(Treatname,Treatsamplenum)))
#group<- factor(c(rep('B',5),rep('C',5),rep('M',5)))
y<- DGEList(counts = data,group = group)
y<- calcNormFactors(y, method = "TMM")
# 查看样本信息
y$samples
norm_counts<- cpm(y, normalized.lib.sizes = TRUE)
# 查看标准化数据
head(norm_counts)
write.csv(data,file = sprintf('%s_data_counts.csv',projectname))
write.csv(norm_counts,file = sprintf('%s_norm_counts.csv',projectname))

####差异基因分析####
colnames(data)
#首先根据分组信息构建试验设计矩阵，分组信息中一定要是对照组在前，处理组在后
dd_CM<- data
colnames(dd_CM)
group_compare<-  factor(c(rep(Controlname,Controlsamplenum),rep(Treatname,Treatsamplenum)),levels = c(Controlname,Treatname))
dgelist <- DGEList(counts =dd_CM, group = group_compare)
#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
#（3）标准化，以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#差异表达基因分析
design <- model.matrix(~group_compare)
#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
write.table(lrt, sprintf('%s_%svs%s_DEGS.txt',projectname,Treatname,Controlname), sep = '\t', col.names = NA, quote = FALSE)

# 筛选差异基因
DEGS<- read.table(sprintf('%s_%svs%s_DEGS.txt',projectname,Treatname,Controlname))
colnames(DEGS)
DEGS$gene<- rownames(DEGS)

DEGS$change<- 'Stable'
# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
if (P=='FDR'){
  cut_off_P = as.numeric(args[12])
  DEGS$change = ifelse(DEGS$FDR< cut_off_P & abs(DEGS$logFC) > cut_off_logFC, 
                       ifelse(DEGS$logFC> cut_off_logFC ,'Up','Down'),
                       'Stable')
}else if (P=='Pvalue'){
  cut_off_P= as.numeric(args[12])
  DEGS$change = ifelse(DEGS$PValue< cut_off_P & abs(DEGS$logFC) > cut_off_logFC, 
                       ifelse(DEGS$logFC> cut_off_logFC ,'Up','Down'),
                       'Stable')
}

head(DEGS)
table(DEGS$change)
write.csv(file = sprintf('%s_%svs%s_change.txt',projectname,Treatname,Controlname),DEGS)

####enrichment####
####UP####
DEGS_UP <- DEGS[DEGS$change=='Up',]
gs = bitr(unique(DEGS_UP$gene), fromType="SYMBOL", toType="ENTREZID", OrgDb=db)
####DEGS-GO_up####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
ego.bp_up <- data.frame(ego.bp@result)
ego.bp_up <- ego.bp_up[ego.bp_up$pvalue<0.05,]
ego.bp_up$Description<- factor(ego.bp_up$Description,levels =ego.bp_up[order(ego.bp_up$Count),]$Description)
write_xlsx(
  data.frame(ego.bp_up),
  path = sprintf("%s-%svs%s_ego.bp_up_P%.4f.xlsx",projectname,Treatname,Controlname,cut_off_P)
)

####DOWN####
DEGS_DOWN <- DEGS[DEGS$change=='Down',]
gs = bitr(unique(DEGS_DOWN$gene), fromType="SYMBOL", toType="ENTREZID", OrgDb=db)
####GO_down####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    

ego.bp_down <- data.frame(ego.bp@result)
ego.bp_down <- ego.bp_down[ego.bp_down$pvalue<0.05,]

ego.bp_down$Description<- factor(ego.bp_down$Description,levels =ego.bp_down[order(ego.bp_down$Count),]$Description)

write_xlsx(
  data.frame(ego.bp_down),
  path = sprintf("%s-%svs%s_ego.bp_down_P%.4f.xlsx",projectname,Treatname,Controlname,cut_off_P)
)











