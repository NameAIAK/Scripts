# 清除工作空间中的所有对象
rm(list=ls())

# 设置工作路径
# setwd("")

# 加载所需的包
library(reshape2)
library(ggplot2)
library(RColorBrewer)
# install.packages('egg')
library(ggalluvial)
library(ggh4x)
library(dplyr)
library(egg)
library(tidyr)

library(readxl) 
# 加载种水平物种丰度表，并设置列名和分隔符
data_raw <- read.csv(
  file = "./83fastp_summary.txt",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t",
  na.strings = "",
  fill=TRUE
)
rownames(data_raw) <- data_raw$samplename

group <- read.table(
  file = "./83group.txt",
  comment.char = "",
  row.names = 1,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t"
)

group <- group %>% arrange(Category1)
data <- data_raw[c(rownames(group)),]
data$group <- group$Category1
library(ggplot2)
data$samplename<- factor(data$samplename,levels =data$samplename)
p_t <- ggplot(data=data,mapping=aes(x=samplename,y=after_filter_readsT,fill=group,group=factor(group)))+
  geom_bar(stat="identity")

p_30 <- ggplot(data=data,mapping=aes(x=samplename,y=q30_rate,fill=group,group=factor(group)))+
  geom_bar(stat="identity")

p_len <- ggplot(data=data,mapping=aes(x=samplename,y=read1_mean_length,fill=group,group=factor(group)))+
  geom_bar(stat="identity")





