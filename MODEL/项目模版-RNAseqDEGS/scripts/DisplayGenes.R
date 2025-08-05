#!/usr/bin/Rscript
####load packages#####
library(ggplot2)
library(cowplot)
library(ggpubr)
#cmd:Rscript DisplayGenes.R GSE188774 test/GeneUP-unique.txt test/GSE188774_norm_counts.csv Control Treat 3 1 3 4 4 5
# args[1] 项目名称：GSE188774
# args[2] 数据库名称缩写：test/GeneUP-unique.txt
# args[3] count文件 test/GSE188774_norm_counts.csv
# args[4] 对照组名称：Control
# args[5] 治疗组/实验组名称：Treat
# args[6] Control样本数：3
# args[7] Control开始数据列：1
# args[8] Treat样本数：3
# args[9] Treat开始数据列：4
# args[10] 图的行数nrow：4
# args[11]  图的列数ncol：5

# 查看操作路径
getwd()
# 获取参数
args <- commandArgs(trailingOnly = TRUE)
projectname <- args[1]#GSE188774
#读取处理后的Genes文本文件
Genesfile <- args[2] # 'test/GeneUP-unique.txt'
Genesdf <- read.csv(Genesfile,col.names = 'genename',header = FALSE)
length(Genesdf$genename)
Geneslist <- Genesdf$genename

#读取标准化数据
norm_countsfile <- args[3] #'test/GSE188774_norm_counts.csv'
norm_counts <- read.csv(norm_countsfile,row.names = 1)
# 对照组名称
Controlname <-args[4] #'Control'
# 治疗组或实验组名称
Treatname <- args[5] #'Treat'
# 对照组样本数量
Controlsamplenum <- as.numeric(args[6])
Controlsamplestartcol <- as.numeric(args[7])#2
Controlsampleendcol <- as.numeric(Controlsamplestartcol)+as.numeric(Controlsamplenum)-1
# 治疗组样本数量
Treatsamplenum <- as.numeric(args[8]) #3
Treatsamplestartcol <- as.numeric(args[9]) #4
Treatsampleendcol <- as.numeric(Treatsamplestartcol)+as.numeric(Treatsamplenum)-1
norm_counts <- norm_counts[,c(Controlsamplestartcol:Controlsampleendcol,Treatsamplestartcol:Treatsampleendcol)]

#准备画图数据
dat_t <- as.data.frame(t(norm_counts))

group <- factor(c(rep(Controlname,Controlsamplenum),rep(Treatname,Treatsamplenum)))
dat_t$group <- group
dat_t$group <- factor(dat_t$group,levels = c(Controlname,Treatname))
comparison <- list(c(Controlname,Treatname))
length(Geneslist)

#作图
p_list <- list()
for (gene in Geneslist){
  print(gene)
  p <- ggplot(dat_t,aes_string(x='group',y=gene))+geom_bar(stat = 'summary',fun=mean,width = 0.5)+stat_summary(fun.data = 'mean_sd',geom = 'errorbar',color='black',width=0.25,size=1,position = position_dodge(.9))+theme_classic()+theme(axis.text = element_text(size = 20))+stat_compare_means(comparisons = comparison,method = 't.test',label = 'p.signif',bracket.size = 1,size=5)+labs(x='',y='',title = gene)+theme(plot.title=element_text(size=20))
  
  p_list[[gene]] <- p 
  
}

nrow <- as.numeric(args[10])
ncol <- as.numeric(args[11])
grid_plot <- plot_grid(plotlist =p_list , nrow = nrow,ncol = ncol)
ggsave(
  plot = grid_plot,
  sprintf("%s-%svs%s-gene-barplot.tsig-norm.pdf",projectname,Treatname,Controlname),
  height = nrow*6,
  width = ncol*4
)

