######火山图#######
getwd()
rm(list=ls())
# setwd('../Genus_EnhancedVolcano/')
#######修改p/log2fc#######
# 设置p值和log2FoldChange阈值
pvalue <- 0.05
log2fc <- 0
contrast = c("Category1", "AB", "C")
# sample_num=39
# g1_num=4
# g2_num=35
# g1_name='C'
# g2_name='AB'

# 加载DESeq2和EnhancedVolcano包
library("DESeq2")
library("EnhancedVolcano")
library('dplyr')  
library(vegan)
library(ggpubr)
# install.packages('dplyr')
# setwd('./module/')
count_tmp <- read.csv(
  file = "./HYXM_16S.G.count.tsv.all.tmp",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t",
  na.strings = "",
  fill=TRUE
)


rownames(count_tmp) <- count_tmp$taxonomy
# 删除第一列OTU
count_tmp <- count_tmp[,-1]
# 删除最后一列taxonomy
df <- count_tmp %>% select(-taxonomy)

original_colnames <- colnames(df) 
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
colnames(df) <- cleaned_colnames  

#####剔除样本#####
# df <- df %>% select(-V350092088_96)
count_tmp <- df
samples <- colnames(df)

map <- read.table(
  file = "HYXM_16S.red.group.txt",
  comment.char = "",
  row.names = 1,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t"
)

# 从样本信息中提取Category1信息
category <- map["Category1"]

# 移除NA值
category <- na.omit(category)
###剔除group样本信息####
category$rowname=rownames(category)
# category <- category %>% filter(rowname != "V350092088_96")
category <- category %>% select(-rowname)

# 匹配分组表和计数表
df <- df[rownames(category)]

# 在计数数据中加1，以避免在后续步骤中对零取对数
df <- df + 1
category$Category1 <- factor(category$Category1)

# 从计数数据和样本信息创建DESeqDataSet对象
dds <- DESeqDataSetFromMatrix(
  countData = df,
  colData = category,
  design = ~Category1 # 比较的分组方案
)

# 将负二项分布模型拟合到数据并进行差异表达测试
fit <- DESeq(dds, quiet = FALSE)


# contrast = c("Category1", "O", "C")
# 提取差异表达分析的结果，将对比设置为将O组与C组进行比较
df <- results(
  fit,
  contrast = contrast
)

# 将结果转换为数据框架，并添加一个Features列，其中包含行名
df <- data.frame(
  Features = rownames(df),
  df,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 设置要在火山图的y轴上使用的列，padj列包含校正后的p值
pcol <- "padj"


# pvalue <- 0.001
# log2fc <- 2
# 设置在火山图中使用的颜色
colors <- c("#00FF00","#FF0000","#666666")

# 按p值和log2FoldChange绝对值对数据框进行排序
df <- df[
  order(
    df[, pcol],
    abs(df$log2FoldChange),
    decreasing = c(FALSE, TRUE)
  ),
]

# 创建一个新列，其中包含p值的-log10
df$`-log10(pvalue)` <- -log(df$pvalue, 10)

# 创建一个新列，其中包含校正后p值的-log10
df$`-log10(padj)` <- -log(df$padj, 10)

# 根据基因的显著性创建点在火山图中的颜色向量
keyvals <- ifelse(
  df$log2FoldChange < -log2fc & df[, pcol] < pvalue,
  colors[1],
  ifelse(
    df$log2FoldChange > log2fc & df[, pcol] < pvalue,
    colors[2],
    colors[3]
  )
)

# 将NA值的颜色设置为灰色
keyvals[is.na(keyvals)] <- colors[3]

# 为图例设置颜色的名称
names(keyvals)[keyvals == colors[1]] <- "Down"
names(keyvals)[keyvals == colors[2]] <- "Up"
names(keyvals)[keyvals == colors[3]] <- "NS"

# 提取前20个基因的标签
# labels <- df[, "Features"]
# select_labs <- head(labels, 20)

# 使用EnhancedVolcano包创建火山图
p_volcano <- EnhancedVolcano(
  df,
  lab = '', # 基因的标签
  xlim=c(-2,2),
  # ylim = c(0,10),
  x = "log2FoldChange", # 用于x轴的变量
  y = pcol, # 用于y轴的变量
  # selectLab = select_labs, # 要突出显示的前几个基因的标签
  xlab = "Log2(Fold Change)", # x轴的标签
  ylab = "-Log10(padj)", # y轴的标签
  pCutoff = pvalue, # 显著性的p值阈值
  FCcutoff = log2fc, # 显著性的log2倍增阈值
  pointSize = 2, # 点的大小
  labSize = 3, # 基因标签的大小
  labCol = "black", # 基因标签的颜色
  title = "AB v.s. C", 
  titleLabSize = 12,# 图表的标题
  subtitle = NULL, # 图表的子标题
  caption = NULL, # 图表的注释
  gridlines.minor = FALSE, # 是否显示次要网格线
  gridlines.major = FALSE, # 是否显示主要网格线
  # boxedLabels = TRUE, # 是否将基因标签放在方框中
  colAlpha = 0.3, # 颜色的alpha值
  legendPosition = "right", # 图例的位置
  legendLabSize = 14, # 图例标签的大小
  legendIconSize = 4.0, # 图例图标的大小
  drawConnectors = TRUE, # 是否绘制连接点和标签的线
  widthConnectors = 1.0, # 连接线的宽度
  colConnectors = "black", # 连接线的颜色
  # shape = c(4, 6, 2, 19), # 点的形状
  # shape = c(19, 19, 19, 19),
  colCustom = keyvals, # 点的颜色
  # 点形状的图例标签
  # legendLabels = c(
  #   sprintf("%s >= %s   and   |log2FC| < %s", pcol, pvalue, log2fc),
  #   sprintf("|log2FC| >= %s  but   %s >= %s", log2fc, pcol, pvalue),
  #   sprintf("%s < %s    but   |log2FC| < %s", pcol, pvalue, log2fc),
  #   sprintf("%s < %s    and   |log2FC| >= %s", pcol, pvalue, log2fc)
  # )
) + coord_cartesian(ylim = c(0,5))+ scale_y_continuous(expand = c(0, 0, 0, 0))+theme(plot.title = element_text(hjust = 0.5)) # 设置图例标签 # 调整y轴的范围


p_volcano

# 将火山图另存为SVG文件
ggsave(
  plot = p_volcano,
  "RHYXM_16S.all.p_volcano.pdf",
  height = 6,
  width = 6
)

######小提琴图########
# p_list <- list() 
# for (i in 1:length(colnames(subset_df)[1:10])) {  
#   variable <- colnames(subset_df)[1:10][i]
#   print(variable) 
#   p <- ggplot(subset_df, aes(x = group, y = .data[[variable]],fill=group)) +  
#     geom_violin(scale = "width", adjust = 0.5) +  # 调整小提琴图的宽度和形状  
#     theme_bw() +  # 使用简洁的主题  
#     labs(title = variable, x = "", y = "abundance") +scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +# 添加标题和轴标签
#     ggtitle(variable)  
#   p_list[[i]] <- p  
# }
# 
# plot_grid(plotlist =p_list , nrow = 2)


library(ggplot2)  

# 创建一个简单的数据集  
set.seed(123)  # 设置随机种子以便结果可复现  
data <- count_tmp
# original_colnames <- colnames(data) 
# cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# # 将修改后的列名赋值回数据框  
# colnames(data) <- cleaned_colnames  
# # names(data)
# 
# data <- data[,-ncol(data)]
group_vlo <- map
group_vlo$sample <- rownames(map)
group_vlo <- group_vlo[,c('sample','Category1')]
colnames(group_vlo) <- c('sample','group')
# group_vlo <- group_vlo %>% filter(sample != "V350092088_96")
# rownames(data) <- data$clade_name
# data
# data <- data[,-1]
data$sum <- rowSums(data)

data1 <- data[order(data$sum, decreasing=TRUE), ]
data1 <- data1[,-ncol(data1)]
# data1 <- data1[1:10,]
data1T <- t(data1)
class(data1T)
data2 <- as.data.frame(data1T)
data2 <- scale(
  data2, 
  center = TRUE, # 减去均值
  scale = TRUE # 除以标准差
)

data2 <- as.data.frame(data2)
data2$group <- group_vlo$group
Genus <- rownames(data1)

# install.packages('cowplot')
library(cowplot)

df_up <- df[df$log2FoldChange > log2fc & df[, pcol] < pvalue,]
df_up <- df_up[complete.cases(df_up), ] 
df_upf <- df_up$Features

df_up <- df_up %>%  
  mutate(change = "UP")  
write.table(df_up, file = "RHYXM_16S.Genus.up.csv", sep = ",", row.names = FALSE, quote = FALSE, na = "")


# df_upf

df_down <- df[df$log2FoldChange < -log2fc & df[, pcol] < pvalue,]
df_down <- df_down[complete.cases(df_down), ] 
df_downf <- df_down$Features
df_down <- df_down %>%  
  mutate(change = "DOWN")
write.table(df_down, file = "RHYXM_16S.Genus.down.csv", sep = ",", row.names = FALSE, quote = FALSE, na = "")

change <- bind_rows(df_up,df_down)
write.table(change, file = "RHYXM_16S.Genus.change.csv", sep = ",", row.names = FALSE, quote = FALSE, na = "")

# df_downf

rbind_df <- rbind(df_up, df_down) 

variables <- rbind_df$Features
length(variables)
# variables <- c('Eubacterium; s__limosum','Lactobacillus; s__johnsonii','Escherichia; s__coli','Shigella; s__dysenteriae','Lactobacillus; s__sp. PV037','Escherichia; s__fergusonii','Shigella; s__sonnei','Escherichia; s__marmotae')
# 使用write函数写入文件（注意：write函数会在每个元素后添加一个换行符）  
write(variables, file = "RHYXM_16S.all.variables.txt") 
p_value_up <- df_up$pvalue[1:length(df_upf)]
p_value_down <- df_down$pvalue[1:length(df_downf)]

# 根据p值确定显著性标记（例如，使用`signif`函数）  
signif_level_up <- ifelse(p_value_up < 0.001, "***", ifelse(p_value_up < 0.01, "**", ifelse(p_value_up < 0.05, "*", "")))
signif_level_down <- ifelse(p_value_down < 0.001, "***", ifelse(p_value_down < 0.01, "**", ifelse(p_value_down < 0.05, "*", "")))

######UP/DOWN不作小提琴图#####

#####单个小提琴图####
# value_var <- c(data2[['Diegoall; s__Yersinia phage phiR1-37']])
# # print(value_var)
# max <- max(value_var) * 0.95
# 
# p <- ggplot(data2, aes(x = group, y = .data[['Diegoall; s__Yersinia phage phiR1-37']],fill=group)) +
#   geom_violin(scale = "width", adjust = 0.5) +  # 调整小提琴图的宽度和形状
#   theme_bw() +   theme(panel.grid.major = element_blank(), # 删除主要网格线  
# panel.grid.minor = element_blank())+# 使用简洁的主题
#   labs(title = 'Diegoall; s__Yersinia phage phiR1-37', x = "", y = "reads") +scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +# 添加标题和轴标签
#   ggtitle('Diegoall; s__Yersinia phage phiR1-37')  +
#   annotate("text", x = 1.5, y = max, label = '', size = 6,
#            vjust = 0, color = "black")
# p

####上调小提琴图####
p_list <- list()
i=1
for (variable in df_upf[1:17]) { 
  
  # print(variable)
  value_var <- c(data2[[variable]])
  # print(value_var)
  max <- max(value_var) * 0.6
  # print(max)
  # class(data2[variable])
  p <- ggplot(data2,aes(x = group, y =.data[[variable]],fill=group))+
    scale_fill_manual(values = c('#00458B','#EE0000'))+
    theme_classic()+
    theme(axis.text = element_text(color = 'black'),legend.position = 'none')+
    geom_violin(trim = F,draw_quantiles = c(0.25, 0.5, 0.75),aes(fill=group) )+
    geom_boxplot(width = 0.1) +
    labs(title = variable, x = "", y = "")+
    geom_point()+
    annotate("text", x = 1.5, y = max, label = signif_level_up[i], size = 6,   
             vjust = 0, color = "black")  # 调整x和y的值以适应你的图
  
  p_list[[variable]] <- p 
  i=i+1
}
 

grid_plot_up <- plot_grid(plotlist =p_list , nrow = 4,ncol = 5)
grid_plot_up
ggsave(
  plot = grid_plot_up,
  "HYXM_16S.all.grid_plot.up.pdf",
  height = 12,
  width = 16
)

# signif_level_up[1:56]
####下调小提琴图####
p_list <- list()
i=1
for (variable in df_downf[1:50]) { 
  
  # print(variable)
  value_var <- c(data2[[variable]])
  # print(value_var)
  max <- max(value_var) * 0.6
  # print(max)
  # class(data2[variable])
  p <- ggplot(data2,aes(x = group, y =.data[[variable]],fill=group))+
    scale_fill_manual(values = c('#00458B','#EE0000'))+
    theme_classic()+
    theme(axis.text = element_text(color = 'black'),legend.position = 'none')+
    geom_violin(trim = F,draw_quantiles = c(0.25, 0.5, 0.75),aes(fill=group) )+
    geom_boxplot(width = 0.1) +
    labs(title = variable, x = "", y = "")+
    geom_point()+
    annotate("text", x = 1.5, y = max, label = signif_level_down[i], size = 6,   
             vjust = 0, color = "black")  # 调整x和y的值以适应你的图
  
  p_list[[variable]] <- p 
  i=i+1
}


grid_plot_down <- plot_grid(plotlist =p_list , nrow = 7,ncol =8)
# grid_plot
ggsave(
  plot = grid_plot_down,
  "HYXM_16S.all.grid_plot.down.pdf",
  height = 12,
  width = 16
)

####合并上调下调小提琴图####
p_list_all <- list()
p_list_all[[1]] <- grid_plot_up
p_list_all[[2]] <- grid_plot_down

grid_plot <- plot_grid(plotlist=p_list_all,ncol = 2,rel_widths = c(6, 1))
# grid_plot
ggsave(
  plot = grid_plot,
  "HYXM_16S.all.grid_plot.all.pdf",
  height = 16,
  width = 32
)

####特别物种####
df_s_up <- df[df$log2FoldChange>0,]
df_s_down <- df[df$log2FoldChange<0,]


index = ( grepl("g__Veillonella", df_s_up$Features)|grepl("g__Barnesiella", df_s_up$Features)|grepl("g__Candidatus_Stoquefichus", df_s_up$Features)|grepl("g__Megamonas", df_s_up$Features))
# df_a = df[index,]
# UP
df_a=df_s_up[index,'Features']
df_a
p_value_a=df_s_up[index,'pvalue']
# DOWN
index = ( grepl("g__Veillonella", df_s_down$Features)|grepl("g__Barnesiella", df_s_down$Features)|grepl("g__Candidatus_Stoquefichus", df_s_down$Features)|grepl("g__Megamonas", df_s_down$Features))

df_a=df_s_down[index,'Features']
df_a
p_value_a=df_s_down[index,'pvalue']

signif_level_a<- ifelse(p_value_a < 0.001, "***", ifelse(p_value_a < 0.01, "**", ifelse(p_value_a < 0.05, "*", "")))
signif_level_a
features <- df_a
signif_level_a <- signif_level_a
p_list <- list()
i=1
for (variable in features) {
  # print(variable)
  # title <- sapply(strsplit(variable, "; g__"), function(x) x[length(x)])
  # title <- gsub("s__", "", title)
  title <- variable
  value_var <- c(data2[[variable]])
  # print(value_var)
  max <- max(value_var) * 0.8
  # print(max)
  # class(data2[variable])
  p <- ggplot(data2,aes(x = group, y =.data[[variable]],fill=group))+
    scale_fill_manual(values = c('#00458B','#EE0000'))+
    theme_classic()+
    theme(axis.text = element_text(color = 'black'),legend.position = 'none')+
    geom_violin(trim = F,draw_quantiles = c(0.25, 0.5, 0.75),aes(fill=group) )+
    geom_boxplot(width = 0.1) +
    labs(title = title, x = "", y = "")+
    geom_point()+
    annotate("text", x = 1.5, y = max, label = signif_level_a[i], size = 6,
             vjust = 0, color = "black")  # 调整x和y的值以适应你的图

  p_list[[variable]] <- p
  i=i+1
  # print(i)
  # signif_level_down[i]
}

grid_plot <- plot_grid(plotlist =p_list , nrow = 2,ncol =3)
grid_plot
ggsave(
  plot = grid_plot,
  "HYXM_16S.cor.down.pdf",
  height = 6,
  width = 8
)
dev.off()

####特别物种箱线图####
# 设置每个组的颜色
group_colors <- c("C" = "#BD3C29", "AB" = "#78D3AC")
# 首先设置比较的列表
compare_list <- list(
  c("C","AB"))
p_s <- ggplot(data2,aes(x = group, y =.data[[variable]],fill=group, colour = group)) +
  geom_boxplot(width = 0.5, alpha = 0.6, lwd = 1.15, outlier.shape = NA) +  # 调整箱的大小
  geom_jitter(width = 0.3, size = 3, alpha = 0.75) +  # 添加散点
  labs(y = "abundance") +
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

p_s
#####热图####
# 安装并加载所需的R包
# install.packages("pheatmap")
library("pheatmap")
# 读入样本分组数据
# map <- read.table(
#   "group.txt",  # 文件名
#   quote = "",  # 引号
#   row.names = 1,  # 第一列作为行名
#   na.strings = "",  # 缺失值表示方式
#   comment.char = "",  # 注释符号
#   check.names = F,  # 是否检查名称
#   stringsAsFactors = F,  # 是否将字符串转换为因子
#   header = TRUE,  # 是否包含表头
#   sep = "\t"  # 分隔符
# )
# 
# # 提取Category1列并去除缺失值
# group <- map["Category1"]
# group <- na.omit(group)
group <- category
# 按行名和Category1列排序
group <- group[order(rownames(group)), , drop = FALSE]
group <- group[order(group[, 1]), , drop = FALSE]

# 读入丰度数据
# abundance <- read.table(
#   "HYXM_16S.S.count.tsv.all.tmp",  # 文件名
#   row.names = 1,  # 第一列作为行名
#   quote = "",  # 引号
#   comment.char = "",  # 注释符号
#   check.names = F,  # 是否检查名称
#   stringsAsFactors = F,  # 是否将字符串转换为因子
#   header = TRUE,  # 是否包含表头
#   sep = "\t"  # 分隔符
# )
# taxonomy <- sapply(strsplit(abundance$taxonomy, "; g__"), function(x) x[length(x)])
# taxonomy <- gsub("; s__", ";", taxonomy)
# rownames(abundance) <- taxonomy

# 去除非数值列并按group列筛选样本
# abundance <- abundance[, sapply(abundance, is.numeric), drop = FALSE]
# original_colnames <- colnames(abundance) 
# cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
# colnames(abundance) <- cleaned_colnames 

# abundance <- abundance[, rownames(group), drop = FALSE]

abundance <- count_tmp
# 转置丰度数据
abundance <- t(abundance)

# 计算各列的总和并选择前20名
# top_idx <- head(
#   order(
#     colSums(abundance),  # 计算各微生物的总丰度
#     decreasing = T  # 按照微生物的总丰度从大到小的顺序排列
#   ),
#   20  # 选择前20名
# )
# 提取前20名的微生物丰度数据
# abundance <- abundance[, c(df_upf[1:20],df_downf[1:20]), drop = FALSE]
abundance <- abundance[, variables, drop = FALSE]

abundance <- abundance + 1
df_log <- log(abundance)
# 标准化（Z-Score）丰度数据的列进行标准化
df_log <- scale(
  df_log,
)
df_log <- df_log[c(rownames(group)),]
# # 计算分组均值
# abundance <- apply(abundance, 2, function(x) {
#   tapply(x, INDEX = group[, 1], mean)
# })


# 定义样本分组颜色
groups_color <- c("#20B2AA", "#DC143C")
names(groups_color) <- unique(group$Category1)
df_log <- t(df_log)

taxonomy <- sapply(strsplit(rownames(df_log), ";g__"), function(x) x[length(x)])
rownames(df_log) <- taxonomy

# # 定义一个逻辑向量，用于指定在哪里插入间隙  
# gap_vec <- rep(FALSE, nrow(mat))  
# # 在第5行和第10行之后插入间隙  
# gap_vec[c(6,9)] <- TRUE  
# 绘制热图
pheatmap(
  df_log,  # 丰度数据
  # scale = 'row',
  legend_breaks = -1:1,
  # file = "HYXM_16S.all.heatmap_cluster.pdf",  # 输出文件名
  file = "HYXM_16S.all.heatmap.pdf",
  fontsize = 10,  # 字体大小
  border_color = "white",  # 边框颜色
  # color = colorRampPalette(  # 颜色渐变
    # colors = c("#20A4B2", "#F9F7F7", "#FF0000")
  # )(100),
  # cluster_cols =TRUE,  # 是否对列进行聚类
  cluster_cols =FALSE,
  clustering_distance_cols = "euclidean",  # 列聚 类距离计算方法
  # cluster_rows = TRUE,  # 是否对行进行聚类
  cluster_rows = FALSE,
  # scale = 'column',
  cellwidth = 10,  # 单元格宽度
  cellheight = 10,  # 单元格高度
  show_rownames = TRUE,  # 是否显示行名
  show_colnames = TRUE,  # 是否显示列名
  annotation_col  = group,  # 行注释
  annotation_colors = list(  # 注释颜色
    Category1=groups_color
  ),
  gaps_row = c(length(df_upf))
)
dev.off()

