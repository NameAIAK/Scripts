#!/usr/bin/Rscript
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

# cmd:Rscript MergeRNAseqDData.R GSE188774 path str1 str2
# args[1] 项目名称：GSE188774
# args[2] 所有count文件夹
# args[3] 样本名称需要删除的字符串
# args[4] 样本名称需要删除的字符串


# 查看操作路径
getwd()
anapath <- getwd()
# 获取参数
args <- commandArgs(trailingOnly = TRUE)
projectname <- args[1]#GSE188774

####合并所有样本数据####
dirpath=args[2] #'dirpath'
setwd(dirpath)
dir=dir()
dat_all <- read.csv(dir[1],sep = '\t')
# head(dat_all)
for (dat in dir[2:length(dir)]){
  print(dat)
  dat_content <- read.csv(dat,sep = '\t')
  dat_all <- merge(dat_all,dat_content,by='Geneid')
}
colnames(dat_all)
replace1 <- args[3]  #"X.data.snakemake.RNA.seq_up_mouse.bowtie2."
replace2 <- args[4] #".bam"
new_colnames <- gsub(replace1, "", names(dat_all))
colnames(dat_all) <- new_colnames
new_colnames <- gsub(replace2, "", names(dat_all))
# 更新数据框的列名
colnames(dat_all) <- new_colnames
setwd(anapath)
write.csv(dat_all,sprintf('%s_dat_all.csv',projectname),row.names = FALSE)


