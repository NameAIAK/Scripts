# install.packages("ggplot2")  
# 假设的数据集  
#设置工作目录
# setwd("")
rm(list=ls())
#####获取alpha_diversity#####
# install.packages('vegan')
library(vegan)
library(ggplot2)
library(dplyr)
otu_table  <-  read.delim("./HYXM_16S.G.count.tsv.all.tmp", header=T, sep="\t", row.names=1, stringsAsFactors = FALSE)
otu_table <- otu_table %>% select(-taxonomy)
original_colnames <- colnames(otu_table) 
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
colnames(otu_table) <- cleaned_colnames  
Transp_otu <- t(otu_table )

Alpha_diversity_index <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Obs <-  est[1, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')
  Pielou <- Shannon / log(Obs, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- rbind(est, Shannon, Simpson,
                  Pielou, goods_coverage)
  if (!is.null(tree)) {
    Pd <- pd(x, tree, include.root = FALSE)[1]
    Pd <- t(Pd)
    result <- rbind(result, Pd)
  }
  result <- as.data.frame(t(result))
  return(result)
}
# est <- estimateR(Transp_otu)
alpha_diversity <- Alpha_diversity_index(Transp_otu)
# alpha_diversity[c(1:3), ]


metadata <- read.delim("./HYXM_16S.red.group.txt", header=T, sep="\t", stringsAsFactors = FALSE)
colnames(metadata) <- c('Sample','Group')

alpha_diversity$Sample <- rownames(alpha_diversity)
plot_data <- merge(metadata, alpha_diversity, by = "Sample")

# custom_theme <- theme(panel.background = element_blank(),
#                       panel.grid.major.y = element_line(colour = "black"), 
#                       axis.line.x = element_line(colour = "black"),
#                       axis.line.y = element_line(colour = "black"),
#                       axis.title.x = element_blank()
# )
# 
# p_richness <- ggplot(data = plot_data, aes(x = Group, y = S.obs))+ 
#   geom_jitter() + 
#   geom_boxplot(aes(fill = Group), alpha = 0.5) + 
#   custom_theme
# ggsave(
#   plot = p_richness,
#   "HYXM_16S.p_richness.pdf",
#   height = 8,
#   width = 10
# )
# 
# p_shannon <- ggplot(data = plot_data, aes(x = Group, y = Shannon))+ 
#   geom_jitter() + 
#   geom_boxplot(aes(fill = Group), alpha = 0.5) + 
#   custom_theme
# ggsave(
#   plot = p_shannon,
#   "HYXM_16S.p_shannon.pdf",
#   height = 8,
#   width = 10
# )
# 
# p_simpson <- ggplot(data = plot_data, aes(x = Group, y = Simpson))+ 
#   geom_jitter() + 
#   geom_boxplot(aes(fill = Group), alpha = 0.5) + 
#   custom_theme
# ggsave(
#   plot = p_simpson,
#   "HYXM_16S.p_simpson.pdf",
#   height = 8,
#   width = 10
# )

#####alpha多样性可视化+差异性######
# install.packages("ggpubr")
library(ggpubr)

# 设置每个组的颜色
group_colors <- c("C" = "#BD3C29", "AB" = "#78D3AC")
# 首先设置比较的列表
compare_list <- list(
  c("C","AB"))

####wilcox-Richness箱线图####
p_Richness <- ggplot(plot_data, aes(x = Group, y = S.obs, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.5, alpha = 0.6, lwd = 1.15, outlier.shape = NA) +  # 调整箱的大小
  geom_jitter(width = 0.3, size = 3, alpha = 0.75) +  # 添加散点
  labs(y = "Richness Diversity") +
  scale_fill_manual(values = group_colors) +  # 设置颜色
  scale_color_manual(values = group_colors) +
  theme_bw() +
  labs(x = NULL) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1.1),
        panel.grid.major = element_blank(), # 去除主网格线
        panel.grid.minor = element_blank() # 去除次网格线
  ) +
  stat_compare_means(
    comparisons = compare_list,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE)# 添加检验结果

# 显示图形
#print(p)
# 保存图形为PDF
ggsave("RHYXM_16S.wilcox.Richness.all.pdf", p_Richness, height = 5, width = 5)


####wilcox-Shannon箱线图####
p_Shannon <- ggplot(plot_data, aes(x = Group, y = Shannon, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.5, alpha = 0.6, lwd = 1.15, outlier.shape = NA) +  # 调整箱的大小
  geom_jitter(width = 0.3, size = 3, alpha = 0.75) +  # 添加散点
  labs(y = "Shannon Diversity") +
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

# 显示图形
#print(p)
# 保存图形为PDF
ggsave("RHYXM_16S.wilcox.Shannon.all.pdf", p_Shannon, height = 5, width = 5)




####wilcox-Simpson箱线图####
p_Simpson <- ggplot(plot_data, aes(x = Group, y = Simpson, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.5, alpha = 0.6, lwd = 1.15, outlier.shape = NA) +  # 调整箱的大小
  geom_jitter(width = 0.3, size = 3, alpha = 0.75) +  # 添加散点
  labs(y = "Simpson Diversity") +
  scale_fill_manual(values = group_colors) +  # 设置颜色
  scale_color_manual(values = group_colors) +
  theme_bw() +
  labs(x = NULL) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1.1),
        panel.grid.major = element_blank(), # 去除主网格线
        panel.grid.minor = element_blank()
        ) +
  stat_compare_means(
    comparisons = compare_list,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE)# 添加检验结果

# 显示图形
#print(p)

# 保存图形为PDF
ggsave("RHYXM_16S.wilcox.Simpson.all.pdf", p_Simpson, height = 5, width = 5)

#####建t-Simpson箱线图####
p_tSimpson <- ggplot(plot_data, aes(x = Group, y = Simpson, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.5, alpha = 0.6, lwd = 1.15, outlier.shape = NA) +  # 调整箱的大小
  geom_jitter(width = 0.3, size = 3, alpha = 0.75) +  # 添加散点
  labs(y = "Simpson Diversity") +
  scale_fill_manual(values = group_colors) +  # 设置颜色
  scale_color_manual(values = group_colors) +
  theme_bw() +
  labs(x = NULL) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1.1),
        panel.grid.major = element_blank(), # 去除主网格线
        panel.grid.minor = element_blank()
        ) +
  stat_compare_means(
    comparisons = compare_list,
    method = "t.test",
    label = "p.signif",
    hide.ns = TRUE)# 添加检验结果

# 显示图形
#print(p)
# 保存图形为PDF
ggsave("HYXM_16S.t.Simpson.all.pdf", p_tSimpson, height = 5, width = 5)

