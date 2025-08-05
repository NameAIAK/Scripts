# 一次性保存工作空间所有对象
save.image(file = "workspace.RData")
# 删除所有对象
rm(list=ls())
# 设置工作路径
setwd()#setwd('E:/work/')
# 查看工作路径
getwd()
# 加载包
library(ggplot2)
# 安装包
install.packages('ggplot2')
BiocManager::install('')#生信包
remotes::install_github('')#github包
# 读取文件
dat_all <- read.csv('',sep = ',')
library(readxl)
data <- read_excel("data.xlsx", sheet = "Sheet1")  # 按名称指定sheet
library(openxlsx)
data <- read.xlsx("data.xlsx", sheet = 1)
write.xlsx(data,'data.xlsx')
data <- read.table("data.txt", header = TRUE, sep = "\t")
library(readr)
data <- read_table("data.txt")

load("data.RData")  # 直接加载到环境（对象名不变）
data <- readRDS("data.rds")  # 适用于单对象读取
save(data, file = "data.RData")  # 存储数据
saveRDS(data, file = "data.rds")  # 单对象存储，推荐

library(jsonlite)
data <- fromJSON("data.json")



