rm(list=ls())
library(tidyr)
library(dplyr)

data <- read.csv(
  file = "./HYXM_16S.G.count.tsv.all.tmp",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t",
  na.strings = "",
  fill=TRUE
)

# 使用separate()函数分割combined_data列，并生成6列新数据  
df_split <- separate(data, col = taxonomy, into = paste0("col", 1:6), sep = ";")  
rownames(df_split) <- df_split$OTU_ID
df_split <- df_split[,-1]


library(dplyr)  

# 假设df是你的数据框，并且包含了所有需要求和的列  

df_summarised <- df_split %>%  
  group_by(col6) %>%  
  summarise(  
    across(starts_with("Black"), ~ sum(.x, na.rm = TRUE), .names = "{col}")  
  )
 
acido_actino_sum <- df_summarised[df_summarised$col2 %in% c(" p__Bacteroidota"," p__Proteobacteria"), -1]   
# 提取p__Aquificota的值  
aquificota_values <- df_summarised[df_summarised$col2 %in% c(" p__Firmicutes"," p__Bacillota"), -1]

# 计算比率  
ratios <- acido_actino_sum / aquificota_values 
rownames(ratios) <- 'ratio'
original_colnames <- colnames(ratios) 
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
colnames(ratios) <- cleaned_colnames
ratios <- t(ratios)
ratios <- as.data.frame(ratios)
class(ratios)
metadata <- read.delim("./group.txt", header=T, sep="\t", stringsAsFactors = FALSE)
colnames(metadata) <- c('Sample','Group')

ratios$Sample <- rownames(ratios)

plot_data <- merge(metadata, ratios, by = "Sample")
library(ggpubr)

# 设置每个组的颜色
group_colors <- c("CONTROL" = "#BD3C29", "OSAHS" = "#78D3AC")
# 首先设置比较的列表
compare_list <- list(
  c("AB","C"))

# 使用ggplot2创建Shannon箱线图
p_ratio <- ggplot(plot_data, aes(x = Group, y = ratio, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.5, alpha = 0.6, lwd = 1.15, outlier.shape = NA) +  # 调整箱的大小
  geom_jitter(width = 0.3, size = 3, alpha = 0.75) +  # 添加散点
  labs(y = "ratio") +
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
p_ratio
# 显示图形
#print(p)
# 保存图形为PDF
ggsave("0613.ratio.all.pdf", p_ratio, height = 5, width = 5)

#####查看物种是否存在####
values1 <- df_summarised[df_summarised$col6 == "g__Megamonas", -1]


rownames(values1) <- c('g__Megamonas')
values1 <- t(values1)
values1 <- as.data.frame(values1)
class(values1)
metadata <- read.delim("./HYXM_16S.black.group.txt", header=T, sep="\t", stringsAsFactors = FALSE)
colnames(metadata) <- c('Sample','Group')

values1$Sample <- rownames(values1)

plot_data <- merge(metadata, values1, by = "Sample")
library(ggpubr)

# 设置每个组的颜色
group_colors <- c("AB" = "#BD3C29", "C" = "#78D3AC")
# 首先设置比较的列表
compare_list <- list(
  c("AB","C"))

# 使用ggplot2创建Shannon箱线图
p_values1 <- ggplot(plot_data, aes(x = Group, y = g__Megamonas, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.5, alpha = 0.4, lwd = 1.15, outlier.shape = NA) +  # 调整箱的大小
  geom_jitter(width = 0.3, size = 3, alpha = 0.75) +  # 添加散点
  labs(y = "g__Megamonas") +
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
p_values1
# 显示图形
#print(p)
# 保存图形为PDF
ggsave("0613.g__Megamonas.all.pdf", p_values1, height = 5, width = 5)
