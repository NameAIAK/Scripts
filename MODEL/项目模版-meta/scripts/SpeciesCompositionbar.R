#下载和加载包
####比对率统计####                       
library("ggplot2") 

count_tmp <- read.csv(
  file = "./osahs.map.stats",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t",
  na.strings = "",
  fill=TRUE
)
 
#绘图
p <- ggplot(count_tmp, aes(x = samplename, y = ratio)) +      
  geom_bar(stat = "identity")+
  theme_classic()+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+ 
  scale_y_continuous(expand=c(0, 0)) + 
  scale_x_discrete(expand=c(0,0))


ggsave(p,filename = 'osahs.map.bar.pdf',height = 8,width = 12)

####病毒占比统计####
rm(list = ls())
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

taxonomy <- sapply(strsplit(data$taxonomy, "; p__"), function(x) x[1])
# rownames(data) <- taxonomy
# 删除第一列OTU
data <- data[,-1]
# 删除最后一列taxonomy
data <- data %>% select(-taxonomy)
df_normalized <- data %>%
  mutate_if(is.numeric, funs(./sum(.))) 
# colSums(df_normalized)
# 添加taxonomy
df_normalized$taxonomy <- taxonomy

original_colnames <- colnames(df_normalized) 
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
colnames(df_normalized) <- cleaned_colnames  

# 计算每一行的和并添加到数据框
library(purrr)  

grouped_df_all <- df_normalized %>%
  group_by(taxonomy) %>%
  summarise_all(sum)

data4 <- melt(df_normalized, id = "taxonomy")
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
p1 <- ggplot(data4, aes(x=variable, y=value, fill = taxonomy)) +
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
p1
p1 <- p1+guides(fill=guide_legend(ncol = 5,byrow = T,
                                  override.aes = list(size = 7),
                                  title.position = "top"))

ggsave(
  plot = p1,
  "Virus.ratio.bar.pdf",
  height = 16,
  width = 20
)








                                         