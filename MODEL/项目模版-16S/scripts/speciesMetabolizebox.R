rm(list=ls())
library(tidyr)
library(dplyr)
library(readxl) 
library(data.table)
####联合分析的样本信息####
num=1
id <- read_excel("HYXM_16S.id.xlsx",sheet = '1')
cor_abs <- 0.3


PMB <- fread(sprintf("../PCancerMb/P-%d.csv",num))
NMB <- fread(sprintf("../PCancerMb/N-%d.csv",num))

id <- na.omit(id)
rownames(id) <- id$meta_id

meta_id <- id[,c("meta_id",'group')]
colnames(meta_id) <- c('sample','group')

P_id <- id[c(meta_id$sample),c("mb_id_P",'group')]
rownames(P_id) <- meta_id$sample
N_id <- id[c(meta_id$sample),c("mb_id_N",'group')]
rownames(N_id) <- meta_id$sample

PMB <- as.data.frame(PMB)
colnames(PMB)[1] <- 'Compound'
# PMB <- PMB[,c(1:50)]
PMB_NAME <- fread(sprintf('../Gred/r0.3/P-%d-差异离子鉴定.csv',num))
# clean_df <- df[!grepl("[^\\x00-\\x7F]", df$text), ]  
df_PMB <- merge(PMB,PMB_NAME,by='Compound') 
# colnames(PMB)
df_PMB <- df_PMB[,c('Accepted Description',P_id$mb_id_P)]

colnames(df_PMB)[1] <- 'name'

df_PMB <- df_PMB  %>%  
  filter(!if_any(c('name'), ~ .x == "" | grepl("^\\s*$", .x)))

mb_pos <- df_PMB[,c('name',P_id$mb_id_P)]
# sort(colnames(PMB))
mb_pos <- mb_pos %>%  
  group_by(name) %>%  
  summarise(  
    across(where(is.numeric), sum, .names = "{.col}"),  
    .groups = 'drop' # 避免后续操作中的分组问题  
  ) 

rownames_mb_pos <- mb_pos$name

mb_pos <- mb_pos[,-1]
rownames(mb_pos) <- rownames_mb_pos
df_mb_pos <- t(mb_pos)


NMB <- as.data.frame(NMB)
colnames(NMB)[1] <- 'Compound'
# NMB <- NMB[,c(1:50)]
NMB_NAME <- fread(sprintf('../Gred/r0.3/N-%d-差异离子鉴定.csv',num))
df_NMB <- merge(NMB,NMB_NAME,by='Compound') 
# colnames(PMB)
df_NMB <- df_NMB[,c('Accepted Description',N_id$mb_id_N)]

colnames(df_NMB)[1] <- 'name'

df_NMB <- df_NMB  %>%  
  filter(!if_any(c('name'), ~ .x == "" | grepl("^\\s*$", .x)))

mb_neg <- df_NMB[,c('name',N_id$mb_id_N)]
# sort(colnames(PMB))
mb_neg <- mb_neg %>%  
  group_by(name) %>%  
  summarise(  
    across(where(is.numeric), sum, .names = "{.col}"),  
    .groups = 'drop' # 避免后续操作中的分组问题  
  ) 


rownames_mb_neg <-mb_neg$name
mb_neg <- mb_neg[,-1]
rownames(mb_neg) <- rownames_mb_neg
df_mb_neg <- t(mb_neg)


df_ts <- df_mb_neg[,c('3-Hydroxykynurenamine')]
df_ts <- as.data.frame(df_ts)
rownames(df_ts) <- rownames(df_mb_neg)
colnames(df_ts) <- c('3-Hydroxykynurenamine')

df_v <- t(df_ts)
values1 <- df_v["3-Hydroxykynurenamine",]
values1 <- as.data.frame(values1)
colnames(values1) <- c('Hydroxykynurenamine')
# values1 <- t(values1)
# values1 <- as.data.frame(values1)
class(values1)
####N/P####
metadata <- N_id
# metadata <- P_id
colnames(metadata)[1] <- 'sample'
# metadata$sample <-rownames(metadata) 
values1$sample <- rownames(values1)

plot_data <- merge(metadata, values1, by = "sample")
library(ggpubr)

# 设置每个组的颜色
group_colors <- c("AB" = "#BD3C29", "C" = "#78D3AC")
# 首先设置比较的列表
compare_list <- list(
  c("AB","C"))

# 使用ggplot2创建Shannon箱线图
p_values1 <- ggplot(plot_data, aes(x = group, y = Hydroxykynurenamine, fill = group, colour = group)) +
  geom_boxplot(width = 0.5, alpha = 0.4, lwd = 1.15, outlier.shape = NA) +  # 调整箱的大小
  geom_jitter(width = 0.3, size = 3, alpha = 0.75) +  # 添加散点
  labs(y = "N7-3-Hydroxykynurenamine") +
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
ggsave("N7-3-Hydroxykynurenamine.pdf", p_values1, height = 5, width = 5)

####读取kraken2分析数据并处理#####
data <- read.csv(
  file = "./HYXM_16S.S.count.tsv.all.tmp",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t",
  na.strings = "",
  fill=TRUE
)
write.csv(data,file ='HYXM_16S.S.count.csv',row.names = FALSE )
count_tmp <- data
count_tmp <- count_tmp[!grepl("g__uncultured", count_tmp$taxonomy), ] 
taxonomy <- sapply(strsplit(count_tmp$taxonomy, ";g__"), function(x) x[length(x)])
taxonomy <- gsub(";s__", ";", taxonomy)
rownames(count_tmp) <- taxonomy
# 删除第一列OTU
count_tmp <- count_tmp[,-1]
# 删除最后一列taxonomy
df <- count_tmp %>% select(-taxonomy)

original_colnames <- colnames(df) 
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
colnames(df) <- cleaned_colnames  


df_abun <- df
samples <- colnames(df)

df_abun <- df_abun[,c(meta_id$sample)]

####SPECIES####
values1 <- df_abun[c('Bifidobacterium;Bifidobacterium_animalis'),]
values1 <- t(values1)
colnames(values1) <- c('Bifidobacterium_Bifidobacterium_animalis')
# values1 <- t(values1)
# values1 <- as.data.frame(values1)
class(values1)

metadata <- meta_id
values1 <- as.data.frame(values1)
values1$sample <- rownames(values1)

plot_data <- merge(metadata, values1, by = "sample")
library(ggpubr)

# 设置每个组的颜色
group_colors <- c("AB" = "#BD3C29", "C" = "#78D3AC")
# 首先设置比较的列表
compare_list <- list(
  c("AB","C"))

# 使用ggplot2创建Shannon箱线图
p_values1 <- ggplot(plot_data, aes(x = group, y = Bifidobacterium_Bifidobacterium_animalis, fill = group, colour = group)) +
  geom_boxplot(width = 0.5, alpha = 0.4, lwd = 1.15, outlier.shape = NA) +  # 调整箱的大小
  geom_jitter(width = 0.3, size = 3, alpha = 0.75) +  # 添加散点
  labs(y = "N7F-Bifidobacterium;Bifidobacterium_animalis") +
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
ggsave("N7F-Bifidobacterium;Bifidobacterium_animalis.pdf", p_values1, height = 5, width = 5)

