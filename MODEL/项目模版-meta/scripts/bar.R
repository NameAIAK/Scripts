# 清除工作空间中的所有对象
rm(list = ls())

# 设置工作路径
# setwd("")

# 加载所需的包
library(reshape2)
library(ggplot2)
library(RColorBrewer)
# install.packages('ggh4x')
library(ggalluvial)
library(ggh4x)
library(dplyr)

# 加载种水平物种丰度表，并设置列名和分隔符
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

taxonomy <- sapply(strsplit(data$taxonomy, "; g__"), function(x) x[length(x)])
rownames(data) <- taxonomy
# 删除第一列OTU
data <- data[,-1]
# 删除最后一列taxonomy
data <- data %>% select(-taxonomy)

original_colnames <- colnames(data) 
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
colnames(data) <- cleaned_colnames  

samples <- colnames(data)


map <- read.table(
  file = "group.txt",
  comment.char = "",
  row.names = 1,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t"
)

# 计算每一行的和并添加到数据框

library(purrr)  
col_sums <- colSums(data) 

df_normalized <- data %>%
  mutate_if(is.numeric, funs(./sum(.))) 
colSums(df_normalized)

# numeric_cols <- sapply(df_normalized, is.numeric) # 找出哪些列是数值型的  
# df_normalized[numeric_cols] <- lapply(df_normalized[numeric_cols], round, digits = 4) 

df_normalized$sum <- rowSums(df_normalized)
# colSums(df_normalized)
# 按每行的和降序排列
data1 <- df_normalized[order(df_normalized$sum, decreasing = TRUE), ]

# 删除sum列
data1 <- data1[, -which(names(data1) == "sum")]
# 取出丰度排名前20的物种，并计算相对丰度
data2 <- data1[1:20, ]

# 计算剩下物种的总丰度
data3 <- 1 - apply(data2, 2, sum)

# 合并数据
data3 <- rbind(data2, data3)
rownames(data3)[nrow(data3)] <- "Others"# 修改最后一行行名others

# 将菌名添加到data3里，为了后面的数据转化
data3$Taxonomy <- factor(rownames(data3), levels = rev(rownames(data3)))

# 宽数据转化为长数据
data4 <- melt(data3, id = "Taxonomy")

group <- map
group$variable <- rownames(group)
names(group)[1] <- 'group'
data4 <- merge(data4, group, by = 'variable')

# 定义调色板
pallet <- c(
  rev(brewer.pal(12, "Paired")), # 反转Paired调色板
  brewer.pal(8, "Set2")[-c(7, 8)], # Set2调色板去除第7和第8个颜色
  brewer.pal(8, "Dark2"), # Dark2调色板
  brewer.pal(12, "Set3"), # Set3调色板
  brewer.pal(8, "Accent"), # Accent调色板
  brewer.pal(11, "Spectral") # Spectral调色板
)
## 无连接线的柱状堆积图
p1 <- ggplot(data4, aes(x=variable, y=value*100, fill = Taxonomy)) +
  #数据输入：样本、物种、丰度
  geom_col(position = 'stack', width = 0.6)+  # stack：堆叠图，width=0.9调整宽度
  scale_y_continuous(expand=c(0, 0))+# 调整y轴属性，使柱子与X轴坐标接触
  scale_fill_manual(values =  pallet) + #手动修改颜色
  
  labs(x = 'Samples', y = 'Relative Abundance(%)') + #设置X轴和Y轴的信息
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        strip.text = element_text(size = 12)) + #设置主题背景，根据自己的需求定制
  # angle：调整横轴标签倾斜角度
  # hjust：上下移动横轴标签
  theme(axis.text.x=element_text(angle=90, hjust=1),axis.text = element_text(size = 12), 
        axis.title = element_text(size = 13), legend.title = element_blank(), 
        legend.text = element_text(size = 11),legend.position="bottom",
        legend.key = element_rect(fill = NA, color = NA)) # 设置边框颜色

p1 <- p1+guides(fill=guide_legend(ncol = 5,byrow = T,
                                  override.aes = list(size = 7),
                                  title.position = "top"))

ggsave(
  plot = p1,
  "0613.all.bar.pdf",
  height = 16,
  width = 20
)

# p1+facet_wrap(~group,scales = 'free_x', ncol = 3)+
#   theme(strip.text = element_text(color = "black", size = 12),# 自定义分面文本
#         strip.background = element_rect(color = "black", fill="grey90"))# 自定义分面背景样式
# 



