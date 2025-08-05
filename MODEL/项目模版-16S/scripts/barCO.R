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

####old####
data_raw <- read.csv(
  file = "./HYXM_16S.G.count.tsv.all.tmp",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t",
  na.strings = "",
  fill=TRUE
)


rownames(data_raw) <- data_raw$taxonomy
# 删除第一列OTU
data <- data_raw[,-1]
# 删除最后一列taxonomy
data <- data %>% select(-taxonomy)

original_colnames <- colnames(data)
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)
# 将修改后的列名赋值回数据框
colnames(data) <- cleaned_colnames

samples <- colnames(data)

map <- read.table(
  file = "./HYXM_16S.red.group.txt",
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
####修改分组信息对应列数####
C_rows <- map %>%  
  filter(Category1 == "C")  

# 如果你需要这些行的行名  
C_rownames <- rownames(C_rows)  

AB_rows <- map %>%  
  filter(Category1 == "AB")  

# 如果你需要这些行的行名  
AB_rownames <- rownames(AB_rows)
data3 <- data3 %>%  
  mutate(C = rowMeans(data3[, c(C_rownames)]),AB=rowMeans(data3[, c(AB_rownames)]))

data3 <- data3[,c('C','AB')]
taxonomy <- sapply(strsplit(rownames(data3), ";g__"), function(x) x[length(x)])
rownames(data3) <- taxonomy

# 将菌名添加到data3里，为了后面的数据转化
data3$Taxonomy <- factor(rownames(data3), levels = rev(rownames(data3)))


# 宽数据转化为长数据
data4 <- melt(data3, id = "Taxonomy")

# group <- map
# group$variable <- rownames(group)
# names(group)[1] <- 'group'
# data4 <- merge(data4, group, by = 'variable')

# 定义调色板
pallet <- c(
  rev(brewer.pal(12, "Paired")), # 反转Paired调色板
  brewer.pal(8, "Set2")[-c(7, 8)], # Set2调色板去除第7和第8个颜色
  brewer.pal(8, "Dark2"), # Dark2调色板
  brewer.pal(12, "Set3"), # Set3调色板
  brewer.pal(8, "Accent"), # Accent调色板
  brewer.pal(11, "Spectral") # Spectral调色板
)

# 使用ggplot2绘制条形图
p <- ggplot(data4, aes(x = variable, y = value)) + # 设置数据和x、y轴变量
  geom_bar(mapping = aes(fill = Taxonomy), stat = "identity",width=0.7) + # 添加条形图层，填充色变量为variable，统计方法为identity，条形宽度为0.7
  theme( legend.key.size = unit(0.2, 'cm'))+theme(legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(title = NULL, ncol = 1, byrow = TRUE),) + # 设置图例标题为空
  scale_fill_manual(values = pallet) + # 设置填充色手动指定调色板
  xlab("") + # 设置x轴标题为空
  ylab("Relative Abundance(%)") + # 设置y轴标题
  theme_bw() + # 设置主题为白色背景
  theme(
    text = element_text(size = 10), # 设置文本大小
    plot.title = element_text(hjust =0.5),
    panel.grid.major = element_blank(), # 去除主网格线
    panel.grid.minor = element_blank(), # 去除次网格线
    axis.line = element_line(), # 设置坐标轴线
    panel.border = element_blank(), # 去除面板边框
    panel.spacing = unit(0, "lines"),      # 面板之间的间距  
    axis.text.x = element_text(angle = 0, size = 6, hjust = 0.5) # 设置x轴文本倾斜角度和大小
  ) +ggtitle('')+
  scale_y_continuous(expand = c(0, 0)) # 设置y轴范围

# 显示条形图
p
ggsave("RHYXM_16S.CO.bar.pdf", egg::set_panel_size( p, width=unit(1, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 300)

data3$diff <- data3$AB-data3$C 
df_sorted <- data3 %>%  
  arrange(desc(diff))  
write.table(df_sorted,file = 'RGenus.top20diff.csv',sep = ',')

rownames_sort <- rownames(df_sorted)
values <- c(df_sorted$diff)
names(values) <- rownames_sort
# barplot(df_sorted$diff, names.arg = df_sorted$Taxonomy, las = 2, main="diff", xlab="", ylab="O-C", col="#1597A5")
# ggsave("HYXM_16S.CO.diff.bar.pdf",p,width = 6, height = 8)
# dev.off()
df_sorted$Taxonomy<- factor(df_sorted$Taxonomy,levels =df_sorted$Taxonomy)

p <- ggplot(df_sorted, aes(x=Taxonomy, y=diff)) +   
  geom_bar(stat="identity", color="black",fill='#1597A5') +  
  theme_minimal() +  
  labs(title="diff", x="", y="AB-C")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  # 将x轴标签旋转90度并调整其大小  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))
ggsave("RHYXM_16S.CO.diff.bar.pdf",p,width = 4, height = 6)

####new####
otu <- data_raw
original_colnames <- colnames(otu) 
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
colnames(otu) <- cleaned_colnames 
otu <- otu %>% select(-taxonomy)
rownames(otu) <- otu$`#OTU ID`
otu <- otu[,-1]

# 使用separate()函数分割combined_data列，并生成6列新数据  
df_split <- separate(data_raw, col = taxonomy, into = paste0("col", 1:7), sep = ";")  
rownames(df_split) <- df_split$`#OTU ID`
df_split <- df_split[,-1]
last_7_cols <- df_split[, -(1:(ncol(df_split) - 7))] 
colnames(last_7_cols) <- c('Kingdom','Phylum',
                           'Class','Order',
                           'Family','Genus',
                           'Species')
tax <- last_7_cols
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")  



library(phyloseq)
# otu_table <-read.csv("otu_table.csv",header = T,row.names = 1)
# 
# taxonomy_table <- read.csv("taxonomy_table.csv",header = T,row.names = 1)
otu_table_phy <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
taxonomy_table_phy <- tax_table(as.matrix(tax))    
#合并
physeq <- phyloseq(otu_table_phy, taxonomy_table_phy) 

####提取Phylum####
phylum_abundance <- tax_glom(physeq, taxrank = "Phylum")
# 提取门的丰度表
phylum_abundance_table <- otu_table(phylum_abundance)
# 提取门的分类信息
phylum_taxonomy <- tax_table(phylum_abundance)
# 将丰度表和分类信息结合
phylum_abundance_with_taxonomy <- cbind(as.data.frame(phylum_abundance_table), as.data.frame(phylum_taxonomy))
# 保存结果到文件
# write.csv(phylum_abundance_with_taxonomy, file = "phylum_abundance_with_taxonomy.csv")

rownames(phylum_abundance_with_taxonomy) <- phylum_abundance_with_taxonomy$Phylum
# 首先获取df的列数  
ncol_df <- ncol(phylum_abundance_with_taxonomy)  
# 然后选择除了最后7列之外的所有列  
phylum_abundance_with_taxonomy <- phylum_abundance_with_taxonomy[, 1:(ncol_df - 7)]
df_abun <- t(phylum_abundance_with_taxonomy)


data <- phylum_abundance_with_taxonomy

map <- read.table(
  file = "./HYXM_16S.group.txt",
  comment.char = "",
  row.names = 1,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t"
)

# 计算每一行的和并添加到数据框
####Phylum物种组成top20####
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
####修改分组信息对应列数####
C_rows <- map %>%  
  filter(Category1 == "C")  

# 如果你需要这些行的行名  
C_rownames <- rownames(C_rows)  

AB_rows <- map %>%  
  filter(Category1 == "AB")  

# 如果你需要这些行的行名  
AB_rownames <- rownames(AB_rows)
data3 <- data3 %>%  
  mutate(C = rowMeans(data3[, c(C_rownames)]),AB=rowMeans(data3[, c(AB_rownames)]))

data3 <- data3[,c('C','AB')]
# 将菌名添加到data3里，为了后面的数据转化
data3$Taxonomy <- factor(rownames(data3), levels = rev(rownames(data3)))

# 宽数据转化为长数据
data4 <- melt(data3, id = "Taxonomy")

# group <- map
# group$variable <- rownames(group)
# names(group)[1] <- 'group'
# data4 <- merge(data4, group, by = 'variable')

# 定义调色板
pallet <- c(
  rev(brewer.pal(12, "Paired")), # 反转Paired调色板
  brewer.pal(8, "Set2")[-c(7, 8)], # Set2调色板去除第7和第8个颜色
  brewer.pal(8, "Dark2"), # Dark2调色板
  brewer.pal(12, "Set3"), # Set3调色板
  brewer.pal(8, "Accent"), # Accent调色板
  brewer.pal(11, "Spectral") # Spectral调色板
)

# 使用ggplot2绘制条形图
p_Phylum <- ggplot(data4, aes(x = variable, y = value)) + # 设置数据和x、y轴变量
  geom_bar(mapping = aes(fill = Taxonomy), stat = "identity",width=0.7) + # 添加条形图层，填充色变量为variable，统计方法为identity，条形宽度为0.7
  theme( legend.key.size = unit(0.2, 'cm'))+theme(legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(title = NULL, ncol = 1, byrow = TRUE),) + # 设置图例标题为空
  scale_fill_manual(values = pallet) + # 设置填充色手动指定调色板
  xlab("") + # 设置x轴标题为空
  ylab("Relative Abundance(%)") + # 设置y轴标题
  theme_bw() + # 设置主题为白色背景
  theme(
    text = element_text(size = 10), # 设置文本大小
    plot.title = element_text(hjust =0.5),
    panel.grid.major = element_blank(), # 去除主网格线
    panel.grid.minor = element_blank(), # 去除次网格线
    axis.line = element_line(), # 设置坐标轴线
    panel.border = element_blank(), # 去除面板边框
    panel.spacing = unit(0, "lines"),      # 面板之间的间距  
    axis.text.x = element_text(angle = 0, size = 16, hjust = 0.5) # 设置x轴文本倾斜角度和大小
  ) +ggtitle('Phylum')+
  scale_y_continuous(expand = c(0, 0)) # 设置y轴范围
  
# 显示条形图
p_Phylum
ggsave("HYXM_16S.Phylum.CO.bar.pdf", egg::set_panel_size( p_Phylum, width=unit(1, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 300)


####提取Class####
Class_abundance <- tax_glom(physeq, taxrank = "Class")
# 提取门的丰度表
Class_abundance_table <- otu_table(Class_abundance)
# 提取门的分类信息
Class_taxonomy <- tax_table(Class_abundance)
# 将丰度表和分类信息结合
Class_abundance_with_taxonomy <- cbind(as.data.frame(Class_abundance_table), as.data.frame(Class_taxonomy))
# 保存结果到文件
# write.csv(Class_abundance_with_taxonomy, file = "Class_abundance_with_taxonomy.csv")

# rownames(Class_abundance_with_taxonomy) <- Class_abundance_with_taxonomy$Class
# 首先获取df的列数  
ncol_df <- ncol(Class_abundance_with_taxonomy)  
# 然后选择除了最后7列之外的所有列  
Class_abundance_with_taxonomy <- Class_abundance_with_taxonomy[, 1:(ncol_df - 7)]
df_abun <- t(Class_abundance_with_taxonomy)


data <- Class_abundance_with_taxonomy

map <- read.table(
  file = "./group.txt",
  comment.char = "",
  row.names = 1,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t"
)

# 计算每一行的和并添加到数据框
####Class物种组成top20####
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
####修改分组信息对应列数####
C_rows <- map %>%  
  filter(Category1 == "C")  

# 如果你需要这些行的行名  
C_rownames <- rownames(C_rows)  

AB_rows <- map %>%  
  filter(Category1 == "AB")  

# 如果你需要这些行的行名  
AB_rownames <- rownames(AB_rows)
data3 <- data3 %>%  
  mutate(C = rowMeans(data3[, c(C_rownames)]),AB=rowMeans(data3[, c(AB_rownames)]))

data3 <- data3[,c('C','AB')]
# 将菌名添加到data3里，为了后面的数据转化
data3$Taxonomy <- factor(rownames(data3), levels = rev(rownames(data3)))

# 宽数据转化为长数据
data4 <- melt(data3, id = "Taxonomy")

# group <- map
# group$variable <- rownames(group)
# names(group)[1] <- 'group'
# data4 <- merge(data4, group, by = 'variable')

# 定义调色板
pallet <- c(
  rev(brewer.pal(12, "Paired")), # 反转Paired调色板
  brewer.pal(8, "Set2")[-c(7, 8)], # Set2调色板去除第7和第8个颜色
  brewer.pal(8, "Dark2"), # Dark2调色板
  brewer.pal(12, "Set3"), # Set3调色板
  brewer.pal(8, "Accent"), # Accent调色板
  brewer.pal(11, "Spectral") # Spectral调色板
)

# 使用ggplot2绘制条形图
p_Class <- ggplot(data4, aes(x = variable, y = value)) + # 设置数据和x、y轴变量
  geom_bar(mapping = aes(fill = Taxonomy), stat = "identity",width=0.7) + # 添加条形图层，填充色变量为variable，统计方法为identity，条形宽度为0.7
  theme( legend.key.size = unit(0.2, 'cm'))+theme(legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(title = NULL, ncol = 1, byrow = TRUE),) + # 设置图例标题为空
  scale_fill_manual(values = pallet) + # 设置填充色手动指定调色板
  xlab("") + # 设置x轴标题为空
  ylab("Relative Abundance(%)") + # 设置y轴标题
  theme_bw() + # 设置主题为白色背景
  theme(
    text = element_text(size = 10), # 设置文本大小
    plot.title = element_text(hjust =0.5),
    panel.grid.major = element_blank(), # 去除主网格线
    panel.grid.minor = element_blank(), # 去除次网格线
    axis.line = element_line(), # 设置坐标轴线
    panel.border = element_blank(), # 去除面板边框
    panel.spacing = unit(0, "lines"),      # 面板之间的间距  
    axis.text.x = element_text(angle = 0, size = 16, hjust = 0.5) # 设置x轴文本倾斜角度和大小
  ) +ggtitle('Class')+
  scale_y_continuous(expand = c(0, 0)) # 设置y轴范围

# 显示条形图
p_Class
ggsave("HYXM_16S.Class.CO.bar.pdf", egg::set_panel_size( p_Class, width=unit(1, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 300)

####提取Order####
Order_abundance <- tax_glom(physeq, taxrank = "Order")
# 提取门的丰度表
Order_abundance_table <- otu_table(Order_abundance)
# 提取门的分类信息
Order_taxonomy <- tax_table(Order_abundance)
# 将丰度表和分类信息结合
Order_abundance_with_taxonomy <- cbind(as.data.frame(Order_abundance_table), as.data.frame(Order_taxonomy))
# 保存结果到文件
# write.csv(Order_abundance_with_taxonomy, file = "Order_abundance_with_taxonomy.csv")

rownames(Order_abundance_with_taxonomy) <- Order_abundance_with_taxonomy$Order
# 首先获取df的列数  
ncol_df <- ncol(Order_abundance_with_taxonomy)  
# 然后选择除了最后7列之外的所有列  
Order_abundance_with_taxonomy <- Order_abundance_with_taxonomy[, 1:(ncol_df - 7)]
df_abun <- t(Order_abundance_with_taxonomy)


data <- Order_abundance_with_taxonomy

map <- read.table(
  file = "./group.txt",
  comment.char = "",
  row.names = 1,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t"
)

# 计算每一行的和并添加到数据框
####Order物种组成top20####
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
####修改分组信息对应列数####
C_rows <- map %>%  
  filter(Category1 == "C")  

# 如果你需要这些行的行名  
C_rownames <- rownames(C_rows)  

AB_rows <- map %>%  
  filter(Category1 == "AB")  

# 如果你需要这些行的行名  
AB_rownames <- rownames(AB_rows)
data3 <- data3 %>%  
  mutate(C = rowMeans(data3[, c(C_rownames)]),AB=rowMeans(data3[, c(AB_rownames)]))

data3 <- data3[,c('C','AB')]
# 将菌名添加到data3里，为了后面的数据转化
data3$Taxonomy <- factor(rownames(data3), levels = rev(rownames(data3)))

# 宽数据转化为长数据
data4 <- melt(data3, id = "Taxonomy")

# group <- map
# group$variable <- rownames(group)
# names(group)[1] <- 'group'
# data4 <- merge(data4, group, by = 'variable')

# 定义调色板
pallet <- c(
  rev(brewer.pal(12, "Paired")), # 反转Paired调色板
  brewer.pal(8, "Set2")[-c(7, 8)], # Set2调色板去除第7和第8个颜色
  brewer.pal(8, "Dark2"), # Dark2调色板
  brewer.pal(12, "Set3"), # Set3调色板
  brewer.pal(8, "Accent"), # Accent调色板
  brewer.pal(11, "Spectral") # Spectral调色板
)

# 使用ggplot2绘制条形图
p_Order <- ggplot(data4, aes(x = variable, y = value)) + # 设置数据和x、y轴变量
  geom_bar(mapping = aes(fill = Taxonomy), stat = "identity",width=0.7) + # 添加条形图层，填充色变量为variable，统计方法为identity，条形宽度为0.7
  theme( legend.key.size = unit(0.2, 'cm'))+theme(legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(title = NULL, ncol = 1, byrow = TRUE),) + # 设置图例标题为空
  scale_fill_manual(values = pallet) + # 设置填充色手动指定调色板
  xlab("") + # 设置x轴标题为空
  ylab("Relative Abundance(%)") + # 设置y轴标题
  theme_bw() + # 设置主题为白色背景
  theme(
    text = element_text(size = 10), # 设置文本大小
    plot.title = element_text(hjust =0.5),
    panel.grid.major = element_blank(), # 去除主网格线
    panel.grid.minor = element_blank(), # 去除次网格线
    axis.line = element_line(), # 设置坐标轴线
    panel.border = element_blank(), # 去除面板边框
    panel.spacing = unit(0, "lines"),      # 面板之间的间距  
    axis.text.x = element_text(angle = 0, size = 16, hjust = 0.5) # 设置x轴文本倾斜角度和大小
  ) +ggtitle('Order')+
  scale_y_continuous(expand = c(0, 0)) # 设置y轴范围

# 显示条形图
p_Order
ggsave("HYXM_16S.Order.CO.bar.pdf", egg::set_panel_size( p_Order, width=unit(1, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 300)

####提取Family####
Family_abundance <- tax_glom(physeq, taxrank = "Family")
# 提取门的丰度表
Family_abundance_table <- otu_table(Family_abundance)
# 提取门的分类信息
Family_taxonomy <- tax_table(Family_abundance)
# 将丰度表和分类信息结合
Family_abundance_with_taxonomy <- cbind(as.data.frame(Family_abundance_table), as.data.frame(Family_taxonomy))
# 保存结果到文件
# write.csv(Family_abundance_with_taxonomy, file = "Family_abundance_with_taxonomy.csv")

rownames(Family_abundance_with_taxonomy) <- Family_abundance_with_taxonomy$Family
# 首先获取df的列数  
ncol_df <- ncol(Family_abundance_with_taxonomy)  
# 然后选择除了最后7列之外的所有列  
Family_abundance_with_taxonomy <- Family_abundance_with_taxonomy[, 1:(ncol_df - 7)]
df_abun <- t(Family_abundance_with_taxonomy)


data <- Family_abundance_with_taxonomy

map <- read.table(
  file = "./group.txt",
  comment.char = "",
  row.names = 1,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t"
)

# 计算每一行的和并添加到数据框
####Family物种组成top20####
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
####修改分组信息对应列数####
C_rows <- map %>%  
  filter(Category1 == "C")  

# 如果你需要这些行的行名  
C_rownames <- rownames(C_rows)  

AB_rows <- map %>%  
  filter(Category1 == "AB")  

# 如果你需要这些行的行名  
AB_rownames <- rownames(AB_rows)
data3 <- data3 %>%  
  mutate(C = rowMeans(data3[, c(C_rownames)]),AB=rowMeans(data3[, c(AB_rownames)]))

data3 <- data3[,c('C','AB')]
# 将菌名添加到data3里，为了后面的数据转化
data3$Taxonomy <- factor(rownames(data3), levels = rev(rownames(data3)))

# 宽数据转化为长数据
data4 <- melt(data3, id = "Taxonomy")

# group <- map
# group$variable <- rownames(group)
# names(group)[1] <- 'group'
# data4 <- merge(data4, group, by = 'variable')

# 定义调色板
pallet <- c(
  rev(brewer.pal(12, "Paired")), # 反转Paired调色板
  brewer.pal(8, "Set2")[-c(7, 8)], # Set2调色板去除第7和第8个颜色
  brewer.pal(8, "Dark2"), # Dark2调色板
  brewer.pal(12, "Set3"), # Set3调色板
  brewer.pal(8, "Accent"), # Accent调色板
  brewer.pal(11, "Spectral") # Spectral调色板
)

# 使用ggplot2绘制条形图
p_Family <- ggplot(data4, aes(x = variable, y = value)) + # 设置数据和x、y轴变量
  geom_bar(mapping = aes(fill = Taxonomy), stat = "identity",width=0.7) + # 添加条形图层，填充色变量为variable，统计方法为identity，条形宽度为0.7
  theme( legend.key.size = unit(0.2, 'cm'))+theme(legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(title = NULL, ncol = 1, byrow = TRUE),) + # 设置图例标题为空
  scale_fill_manual(values = pallet) + # 设置填充色手动指定调色板
  xlab("") + # 设置x轴标题为空
  ylab("Relative Abundance(%)") + # 设置y轴标题
  theme_bw() + # 设置主题为白色背景
  theme(
    text = element_text(size = 10), # 设置文本大小
    plot.title = element_text(hjust =0.5),
    panel.grid.major = element_blank(), # 去除主网格线
    panel.grid.minor = element_blank(), # 去除次网格线
    axis.line = element_line(), # 设置坐标轴线
    panel.border = element_blank(), # 去除面板边框
    panel.spacing = unit(0, "lines"),      # 面板之间的间距  
    axis.text.x = element_text(angle = 0, size = 16, hjust = 0.5) # 设置x轴文本倾斜角度和大小
  ) +ggtitle('Family')+
  scale_y_continuous(expand = c(0, 0)) # 设置y轴范围

# 显示条形图
p_Family
ggsave("HYXM_16S.Family.CO.bar.pdf", egg::set_panel_size( p_Family, width=unit(1, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 300)


####提取Genus####
Genus_abundance <- tax_glom(physeq, taxrank = "Genus")
# 提取门的丰度表
Genus_abundance_table <- otu_table(Genus_abundance)
# 提取门的分类信息
Genus_taxonomy <- tax_table(Genus_abundance)
# 将丰度表和分类信息结合
Genus_abundance_with_taxonomy <- cbind(as.data.frame(Genus_abundance_table), as.data.frame(Genus_taxonomy))
# 保存结果到文件
# write.csv(Genus_abundance_with_taxonomy, file = "Genus_abundance_with_taxonomy.csv")

rownames(Genus_abundance_with_taxonomy) <- Genus_abundance_with_taxonomy$Genus
# 首先获取df的列数  
ncol_df <- ncol(Genus_abundance_with_taxonomy)  
# 然后选择除了最后7列之外的所有列  
Genus_abundance_with_taxonomy <- Genus_abundance_with_taxonomy[, 1:(ncol_df - 7)]
df_abun <- t(Genus_abundance_with_taxonomy)


data <- Genus_abundance_with_taxonomy

map <- read.table(
  file = "./group.txt",
  comment.char = "",
  row.names = 1,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t"
)

# 计算每一行的和并添加到数据框
####Genus物种组成top20####
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
####修改分组信息对应列数####
C_rows <- map %>%  
  filter(Category1 == "C")  

# 如果你需要这些行的行名  
C_rownames <- rownames(C_rows)  

AB_rows <- map %>%  
  filter(Category1 == "AB")  

# 如果你需要这些行的行名  
AB_rownames <- rownames(AB_rows)
data3 <- data3 %>%  
  mutate(C = rowMeans(data3[, c(C_rownames)]),AB=rowMeans(data3[, c(AB_rownames)]))

data3 <- data3[,c('C','AB')]
# 将菌名添加到data3里，为了后面的数据转化
data3$Taxonomy <- factor(rownames(data3), levels = rev(rownames(data3)))

# 宽数据转化为长数据
data4 <- melt(data3, id = "Taxonomy")

# group <- map
# group$variable <- rownames(group)
# names(group)[1] <- 'group'
# data4 <- merge(data4, group, by = 'variable')

# 定义调色板
pallet <- c(
  rev(brewer.pal(12, "Paired")), # 反转Paired调色板
  brewer.pal(8, "Set2")[-c(7, 8)], # Set2调色板去除第7和第8个颜色
  brewer.pal(8, "Dark2"), # Dark2调色板
  brewer.pal(12, "Set3"), # Set3调色板
  brewer.pal(8, "Accent"), # Accent调色板
  brewer.pal(11, "Spectral") # Spectral调色板
)

# 使用ggplot2绘制条形图
p_Genus <- ggplot(data4, aes(x = variable, y = value)) + # 设置数据和x、y轴变量
  geom_bar(mapping = aes(fill = Taxonomy), stat = "identity",width=0.7) + # 添加条形图层，填充色变量为variable，统计方法为identity，条形宽度为0.7
  theme( legend.key.size = unit(0.2, 'cm'))+theme(legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(title = NULL, ncol = 1, byrow = TRUE),) + # 设置图例标题为空
  scale_fill_manual(values = pallet) + # 设置填充色手动指定调色板
  xlab("") + # 设置x轴标题为空
  ylab("Relative Abundance(%)") + # 设置y轴标题
  theme_bw() + # 设置主题为白色背景
  theme(
    text = element_text(size = 10), # 设置文本大小
    plot.title = element_text(hjust =0.5),
    panel.grid.major = element_blank(), # 去除主网格线
    panel.grid.minor = element_blank(), # 去除次网格线
    axis.line = element_line(), # 设置坐标轴线
    panel.border = element_blank(), # 去除面板边框
    panel.spacing = unit(0, "lines"),      # 面板之间的间距  
    axis.text.x = element_text(angle = 0, size = 16, hjust = 0.5) # 设置x轴文本倾斜角度和大小
  ) +ggtitle('Genus')+
  scale_y_continuous(expand = c(0, 0)) # 设置y轴范围

# 显示条形图
p_Genus
ggsave("HYXM_16S.Genus.CO.bar.pdf", egg::set_panel_size( p_Genus, width=unit(1, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 300)

####提取Species####
Species_abundance <- tax_glom(physeq, taxrank = "Species")
# 提取门的丰度表
Species_abundance_table <- otu_table(Species_abundance)
# 提取门的分类信息
Species_taxonomy <- tax_table(Species_abundance)
# 将丰度表和分类信息结合
Species_abundance_with_taxonomy <- cbind(as.data.frame(Species_abundance_table), as.data.frame(Species_taxonomy))
# 保存结果到文件
# write.csv(Species_abundance_with_taxonomy, file = "Species_abundance_with_taxonomy.csv")

rownames(Species_abundance_with_taxonomy) <- Species_abundance_with_taxonomy$Species
# 首先获取df的列数  
ncol_df <- ncol(Species_abundance_with_taxonomy)  
# 然后选择除了最后7列之外的所有列  
Species_abundance_with_taxonomy <- Species_abundance_with_taxonomy[, 1:(ncol_df - 7)]
df_abun <- t(Species_abundance_with_taxonomy)


data <- Species_abundance_with_taxonomy

map <- read.table(
  file = "./group.txt",
  comment.char = "",
  row.names = 1,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t"
)

# 计算每一行的和并添加到数据框
####Species物种组成top20####
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
####修改分组信息对应列数####
C_rows <- map %>%  
  filter(Category1 == "C")  

# 如果你需要这些行的行名  
C_rownames <- rownames(C_rows)  

AB_rows <- map %>%  
  filter(Category1 == "AB")  

# 如果你需要这些行的行名  
AB_rownames <- rownames(AB_rows)
data3 <- data3 %>%  
  mutate(C = rowMeans(data3[, c(C_rownames)]),AB=rowMeans(data3[, c(AB_rownames)]))

data3 <- data3[,c('C','AB')]
# 将菌名添加到data3里，为了后面的数据转化
data3$Taxonomy <- factor(rownames(data3), levels = rev(rownames(data3)))

# 宽数据转化为长数据
data4 <- melt(data3, id = "Taxonomy")

# group <- map
# group$variable <- rownames(group)
# names(group)[1] <- 'group'
# data4 <- merge(data4, group, by = 'variable')

# 定义调色板
pallet <- c(
  rev(brewer.pal(12, "Paired")), # 反转Paired调色板
  brewer.pal(8, "Set2")[-c(7, 8)], # Set2调色板去除第7和第8个颜色
  brewer.pal(8, "Dark2"), # Dark2调色板
  brewer.pal(12, "Set3"), # Set3调色板
  brewer.pal(8, "Accent"), # Accent调色板
  brewer.pal(11, "Spectral") # Spectral调色板
)

# 使用ggplot2绘制条形图
p_Species <- ggplot(data4, aes(x = variable, y = value)) + # 设置数据和x、y轴变量
  geom_bar(mapping = aes(fill = Taxonomy), stat = "identity",width=0.7) + # 添加条形图层，填充色变量为variable，统计方法为identity，条形宽度为0.7
  theme( legend.key.size = unit(0.2, 'cm'))+theme(legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(title = NULL, ncol = 1, byrow = TRUE),) + # 设置图例标题为空
  scale_fill_manual(values = pallet) + # 设置填充色手动指定调色板
  xlab("") + # 设置x轴标题为空
  ylab("Relative Abundance(%)") + # 设置y轴标题
  theme_bw() + # 设置主题为白色背景
  theme(
    text = element_text(size = 10), # 设置文本大小
    plot.title = element_text(hjust =0.5),
    panel.grid.major = element_blank(), # 去除主网格线
    panel.grid.minor = element_blank(), # 去除次网格线
    axis.line = element_line(), # 设置坐标轴线
    panel.border = element_blank(), # 去除面板边框
    panel.spacing = unit(0, "lines"),      # 面板之间的间距  
    axis.text.x = element_text(angle = 0, size = 16, hjust = 0.5) # 设置x轴文本倾斜角度和大小
  ) +ggtitle('Species')+
  scale_y_continuous(expand = c(0, 0)) # 设置y轴范围

# 显示条形图
p_Species
ggsave("HYXM_16S.Species.CO.bar.pdf", egg::set_panel_size( p_Species, width=unit(1, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 300)




# p1+facet_wrap(~group,scales = 'free_x', ncol = 3)+
#   theme(strip.text = element_text(color = "black", size = 12),# 自定义分面文本
#         strip.background = element_rect(color = "black", fill="grey90"))# 自定义分面背景样式
#####合并####
# p_list <- list()
# p_list[['Phylum']] <- p_Phylum
# p_list[['Phylum']] <- p_Phylum
# p_list[['Phylum']] <- p_Phylum
# p_list[['Phylum']] <- p_Phylum
# 
# 
# library(cowplot)
# grid_plot <- plot_grid(plotlist =p_list , nrow = 1,ncol =4)
# # grid_plot
# ggsave(
#   plot = grid_plot,
#   "HYXM_16S.join.pdf",
#   height = 8,
#   width = 40
# )
