######火山图#######
getwd()
rm(list=ls())
# setwd('../KO_volcano/')
#######修改p/log2fc#######
# 设置p值和log2FoldChange阈值
pvalue <- 0.05
log2fc <- 0
####实验组在前####
contrast = c("Category1", "OSAHS", "CONTROL")
# sample_num=39
# g1_num=4
# g2_num=35
# g1_name='control'
# g2_name='osahs'

# 加载DESeq2和EnhancedVolcano包
library("DESeq2")
library("EnhancedVolcano")
library('dplyr')  
# install.packages('dplyr')
# setwd('./module/')
count_tmp <- read.csv(
  file = "./0828.KO_samples.xls",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t",
  na.strings = "",
  fill=TRUE
)

# KO_name <- sapply(strsplit(count_tmp$KO_name, "; g__"), function(x) x[length(x)])
class(count_tmp)
# install.packages("data.table")
library(data.table)  
setDT(count_tmp) # 将df转换为data.table格式  
# 使用unique  
df_unique <- unique(count_tmp, by = "KO", fromLast = FALSE)  
count_tmp <- as.data.frame(df_unique)
count_tmp <- na.omit(count_tmp)
rownames(count_tmp) <- count_tmp$KO

count_tmp <- count_tmp[,-c(1,2,3)]
# df <- count_tmp %>% select(-KO_des)

# original_colnames <- colnames(df) 
# cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
# colnames(df) <- cleaned_colnames  

samples <- colnames(count_tmp)

map <- read.table(
  file = "0828.group.txt",
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

# 匹配分组表和计数表
count_tmp <- count_tmp[rownames(category)]

# 在计数数据中加1，以避免在后续步骤中对零取对数
df <- count_tmp
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
labels <- df[, "Features"]
select_labs <- head(labels, 20)

# 使用EnhancedVolcano包创建火山图
p_volcano <- EnhancedVolcano(
  df,
  lab = '', # 基因的标签
  x = "log2FoldChange", # 用于x轴的变量
  y = pcol, # 用于y轴的变量
  xlim = c(-5,5 ),
  ylim = c(0, 8),
  # selectLab = select_labs, # 要突出显示的前几个基因的标签
  xlab = "Log2(Fold Change)", # x轴的标签
  ylab = "-Log10(padj)", # y轴的标签
  pCutoff = pvalue, # 显著性的p值阈值
  FCcutoff = log2fc, # 显著性的log2倍增阈值
  pointSize = 2, # 点的大小
  labSize = 3, # 基因标签的大小
  labCol = "black", # 基因标签的颜色
  title = "osahs v.s. control", 
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
) +theme(plot.title = element_text(hjust = 0.5))# 设置图例标签 # 调整y轴的范围

p_volcano
# 将火山图另存为SVG文件
ggsave(
  plot = p_volcano,
  "0828.KO2.p_volcano.pdf",
  height = 8,
  width = 10
)

######小提琴图########
# p_list <- list() 
# for (i in 1:length(colnames(subset_df)[1:10])) {  
#   variable <- colnames(subset_df)[1:10][i]
#   print(variable) 
#   p <- ggplot(subset_df, aes(x = group, y = .data[[variable]],fill=group)) +  
#     geom_violin(scale = "width", adjust = 0.5) +  # 调整小提琴图的宽度和形状  
#     theme_minimal() +  # 使用简洁的主题  
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
# 将修改后的列名赋值回数据框  
# colnames(data) <- cleaned_colnames  
# names(data)
# data <- data %>% select(-KO_des)
# group_vlo <- data.frame(Sample = c(samples),group = c(rep(g1_name,g1_num),rep(g2_name,g2_num)))

group_vlo <- map
group_vlo$sample <- rownames(map)
group_vlo <- group_vlo[,c('sample','Category1')]
colnames(group_vlo) <- c('sample','group')
# rownames(data) <- data$clade_name
# data
# data <- data[,-1]
data$sum <- rowSums(data)

data1 <- data[order(data$sum, decreasing=TRUE), ]
data1 <- select(data1,-'sum')
# data1 <- data1[1:10,]
data1T <- t(data1)
class(data1T)
data2 <- as.data.frame(data1T)
data2 <- log(data2)
data2 <- scale(
  data2, 
  center = TRUE, # 减去均值
  scale = TRUE # 除以标准差
)

data2 <- as.data.frame(data2)
data2$group <- group_vlo$group
KOs <- rownames(data1)

# install.packages('cowplot')
library(cowplot)

df_up <- df[df$log2FoldChange > log2fc & df[, pcol] < pvalue,]
df_up <- df_up[complete.cases(df_up), ] 
df_upf <- df_up$Features
df_upf

df_down <- df[df$log2FoldChange < -log2fc & df[, pcol] < pvalue,]
df_down <- df_down[complete.cases(df_down), ] 
df_downf <- df_down$Features
df_downf

rbind_df <- rbind(df_up, df_down) 

variables <- rbind_df$Features
length(variables)
# variables <- c('Eubacterium; s__limosum','Lactobacillus; s__johnsonii','Escherichia; s__coli','Shigella; s__dysenteriae','Lactobacillus; s__sp. PV037','Escherichia; s__fergusonii','Shigella; s__sonnei','Escherichia; s__marmotae')
# 使用write函数写入文件（注意：write函数会在每个元素后添加一个换行符）  
write(variables, file = "0828.KO.variables.txt") 
p_value_up <- df_up$pvalue[1:length(df_upf)]
p_value_down <- df_down$pvalue[1:length(df_downf)]
# 根据p值确定显著性标记（例如，使用`signif`函数）  
signif_level_up <- ifelse(p_value_up < 0.001, "***", ifelse(p_value_up < 0.01, "**", ifelse(p_value_up < 0.05, "*", "")))
signif_level_down <- ifelse(p_value_down < 0.001, "***", ifelse(p_value_down < 0.01, "**", ifelse(p_value_down < 0.05, "*", "")))

# p_list <- list()
# i=1
# for (variable in df_upf) { 
#   # print(variable)
#   value_var <- c(data2[[variable]])
#   # print(value_var)
#   max <- max(value_var) * 0.95
#   # print(max)
#   # class(data2[variable])
#   p <- ggplot(data2, aes(x = group, y = .data[[variable]],fill=group)) +  
#     geom_violin(scale = "width", adjust = 0.5) +  # 调整小提琴图的宽度和形状  
#     theme_minimal() +  # 使用简洁的主题  
#     labs(title = variable, x = "", y = "reads") +scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +# 添加标题和轴标签
#     ggtitle(variable)  +
#     annotate("text", x = 1.5, y = max, label = signif_level_up[i], size = 6,   
#              vjust = 0, color = "black")  # 调整x和y的值以适应你的图
#   p_list[[variable]] <- p 
#   i=i+1
# }
#  
# 
# grid_plot_up <- plot_grid(plotlist =p_list , nrow = 7,ncol = 2)
# # grid_plot_up
# ggsave(
#   plot = grid_plot_up,
#   "0828.KO.grid_plot.up.pdf",
#   height = 50,
#   width = 50,
#   limitsize = FALSE
# )
# ######单个小提琴图####
# # value_var <- c(data2[[variable]])
# # # print(value_var)
# # max <- max(value_var) * 0.95
# # p <- ggplot(data2, aes(x = group, y = .data[['Morganella; s__morganii']],fill=group)) +  
# #   geom_violin(scale = "width", adjust = 0.5) +  # 调整小提琴图的宽度和形状  
# #   theme_minimal() +  # 使用简洁的主题  
# #   labs(title = 'Morganella; s__morganii', x = "", y = "reads") +scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +# 添加标题和轴标签
# #   ggtitle('Morganella; s__morganii')  +
# #   annotate("text", x = 1.5, y = max, label = '*', size = 6,   
# #            vjust = 0, color = "black") 
# # p
# p_list <- list()
# i=1
# for (variable in df_downf) { 
#   # print(variable)
#   value_var <- c(data2[[variable]])
#   # print(value_var)
#   max <- max(value_var) * 0.95
#   # print(max)
#   # class(data2[variable])
#   p <- ggplot(data2, aes(x = group, y = .data[[variable]],fill=group)) +  
#     geom_violin(scale = "width", adjust = 0.5) +  # 调整小提琴图的宽度和形状  
#     theme_minimal() +  # 使用简洁的主题  
#     labs(title = variable, x = "", y = "reads") +scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +# 添加标题和轴标签
#     ggtitle(variable)  +
#     annotate("text", x = 1.5, y = max, label = signif_level_down[i], size = 6,   
#              vjust = 0, color = "black")  # 调整x和y的值以适应你的图
#   p_list[[variable]] <- p 
#   i=i+1
#   # print(i)
#   # signif_level_down[i]
# }
# 
# 
# grid_plot_down <- plot_grid(plotlist =p_list , nrow = 7,ncol =9)
# # grid_plot
# ggsave(
#   plot = grid_plot_down,
#   "0828.KO.grid_plot.down.pdf",
#   height = 8,
#   width = 16
# )
# 
# p_list_all <- list()
# p_list_all[[1]] <- grid_plot_up
# p_list_all[[2]] <- grid_plot_down
# ######修改小提琴all图片大小和占比#######
# grid_plot <- plot_grid(plotlist=p_list_all,ncol = 2,rel_widths = c(2, 9))
# # grid_plot
# ggsave(
#   plot = grid_plot,
#   "0828.KO.grid_plot.all.pdf",
#   height = 30,
#   width = 48
# )
#####热图####
# 安装并加载所需的R包
# install.packages("pheatmap")
library("pheatmap")
# 读入样本分组数据
map <- read.table(
  "0828.group.txt",  # 文件名
  quote = "",  # 引号
  row.names = 1,  # 第一列作为行名
  na.strings = "",  # 缺失值表示方式
  comment.char = "",  # 注释符号
  check.names = F,  # 是否检查名称
  stringsAsFactors = F,  # 是否将字符串转换为因子
  header = TRUE,  # 是否包含表头
  sep = "\t"  # 分隔符
)

# 提取Category1列并去除缺失值
group <- map["Category1"]
group <- na.omit(group)

# 按行名和Category1列排序
group <- group[order(rownames(group)), , drop = FALSE]
group <- group[order(group[, 1]), , drop = FALSE]

# 读入丰度数据
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
# 提取差异物种的微生物丰度数据
class(df_up$Features)
ts_variables=c(df_up$Features[1:10],df_down$Features[1:10])
abundance <- abundance[, ts_variables, drop = FALSE]
abundance <- abundance[, variables, drop = FALSE]
abundance <- abundance + 1
# 标准化（Z-Score）丰度数据-列
abundance <- log(abundance)
abundance <- scale(
  abundance
)
df_log <- abundance[c(rownames(group)),]
# # 计算分组均值
# abundance <- apply(abundance, 2, function(x) {
#   tapply(x, INDEX = group[, 1], mean)
# })


# 定义样本分组颜色
groups_color <- c("#20B2AA", "#DC143C")
names(groups_color) <- unique(group$Category1)
df_log <- t(df_log)

# # 定义一个逻辑向量，用于指定在哪里插入间隙  
# gap_vec <- rep(FALSE, nrow(mat))  
# # 在第5行和第10行之后插入间隙  
# gap_vec[c(6,9)] <- TRUE  
# 绘制热图
pheatmap(
  df_log,  # 丰度数据
  file = "0828.KO2-10.heatmap.pdf",  # 输出文件名
  fontsize = 10,  # 字体大小
  border_color = "white",  # 边框颜色
  # color = colorRampPalette(  # 颜色渐变
  #   colors = c("#20A4B2", "#F9F7F7", "#FF0000")
  # )(100),
  cluster_cols = FALSE,  # 是否对列进行聚类
  clustering_distance_cols = "euclidean",  # 列聚类距离计算方法
  cluster_rows = FALSE,  # 是否对行进行聚类
  # scale = 'row',
  cellwidth = 10,  # 单元格宽度
  cellheight = 10,  # 单元格高度
  show_rownames = TRUE,  # 是否显示行名
  show_colnames = TRUE,  # 是否显示列名
  annotation_col  = group,  # 行注释
  annotation_colors = list(  # 注释颜色
    Category1=groups_color
  ),
######修改热图分割位置#######
  # gaps_row = c(length(df_upf))
  gaps_row = c(10)
)


######KO富集分析######
#enrichment analysis using clusterprofiler package
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
# library(DOSE)
# BiocManager::install("org.Mm.eg.db")
# packagedir <- choose.files('../../org.Mm.eg.db_3.19.1.tar.gz')
# install.packages('../../org.Mm.eg.db_3.19.1.tar.gz' , repos=NULL)
# library(org.Mm.eg.db)
#get the ENTREZID for the next analysis
write.csv(df, "0828.KO_df.csv", row.names = FALSE)
# library("clusterProfiler")
## 将Description转换为因子：
convert_fraction_to_numeric <- function(fraction_str) {  
  # 使用strsplit分割字符串，结果为列表  
  # 使用[[1]]获取列表的第一个元素（这里假设每个字符串只包含一个分数）  
  # 然后用as.numeric转换每个元素为数值  
  parts <- as.numeric(unlist(strsplit(fraction_str, "/")))  
  # 执行除法  
  return(parts[1] / parts[2])  
}
#####导入KOID/UP####
KO_list_up=df_up$Features
#进行富集分析
result_up=enrichKEGG(KO_list_up,
                   organism = "ko", #物种选择ko
                   pvalueCutoff = 0.05,  #p值cutoff
                   pAdjustMethod = "BH", #FDR矫正p值
                   qvalueCutoff = 0.2,   #q值cutoff
)
#保存富集分析结果
write.csv(file="0828.KO_enrichment_result_up.csv",data.frame(result_up@result),row.names=F)

KEGG_result_up <- data.frame(result_up@result)
KEGG_result_up <- KEGG_result_up[KEGG_result_up$pvalue<0.05,]
KEGG_result_up <- KEGG_result_up[c(2:11),]
KEGG_result_up$GeneRatio_numeric <- sapply(KEGG_result_up$GeneRatio, convert_fraction_to_numeric)  
KEGG_result_up[order(KEGG_result_up$GeneRatio_numeric),]$Description
KEGG_result_up$Description<- factor(KEGG_result_up$Description,levels =KEGG_result_up[order(KEGG_result_up$GeneRatio_numeric),]$Description)
write.csv(file="0828.KO_enrichment_result_up_p0.05.csv",data.frame(KEGG_result_up),row.names=F)

# 绘制条形图
pdf(file="0828.KO_enrichment_up_bar.pdf",width = 10,height = 7)

ggplot(data = KEGG_result_up,
            aes(x = GeneRatio_numeric, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "KEGG_UP") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()
# 绘制气泡图
pdf(file="0828.KO_enrichment_up_dot.pdf",width = 10,height = 7)
dotplot(result_up, showCategory = 8)
dev.off()

####导入KOID/DOWN####
KO_list_down=df_down$Features
#进行富集分析
result_down=enrichKEGG(KO_list_down,
                     organism = "ko", #物种选择ko
                     pvalueCutoff = 0.05,  #p值cutoff
                     pAdjustMethod = "BH", #FDR矫正p值
                     qvalueCutoff = 0.2   #q值cutoff
)
#保存富集分析结果
write.csv(file="0828.KO_enrichment_result_down.csv",data.frame(result_down@result),row.names=F)

KEGG_result_down <- data.frame(result_down@result)
KEGG_result_down <- KEGG_result_down[KEGG_result_down$pvalue<0.05,]
KEGG_result_down <- KEGG_result_down[c(1:10),]
KEGG_result_down$GeneRatio_numeric <- sapply(KEGG_result_down$GeneRatio, convert_fraction_to_numeric)  
KEGG_result_down[order(KEGG_result_down$GeneRatio_numeric),]$Description

KEGG_result_down$Description<- factor(KEGG_result_down$Description,levels =KEGG_result_down[order(KEGG_result_down$GeneRatio_numeric),]$Description)
write.csv(file="0828.KO_enrichment_result_down_p0.05.csv",data.frame(KEGG_result_down),row.names=F)

# 绘制条形图
pdf(file="0828.KO_enrichment_down_bar.pdf",width = 10,height = 7)

ggplot(data = KEGG_result_down,
       aes(x = GeneRatio_numeric, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "KEGG_DOWN") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()
# 绘制气泡图
pdf(file="0828.KO_enrichment_down_dot.pdf",width = 10,height = 7)
dotplot(result_down, showCategory = 8)
dev.off()

#####导入KOID/all####
KO_list_all=df$Features
#进行富集分析
result_all=enrichKEGG(KO_list_all,
                       organism = "ko", #物种选择ko
                       pvalueCutoff = 0.05,  #p值cutoff
                       pAdjustMethod = "BH", #FDR矫正p值
                       qvalueCutoff = 0.2   #q值cutoff
)
#保存富集分析结果
write.csv(file="0828.KO_enrichment_result_all.csv",data.frame(result_all),row.names=F)

# 绘制条形图
pdf(file="0828.KO_enrichment_all_bar.pdf",width = 10,height = 7)
class(result_all)
barplot(result_all, x='GeneRatio',drop = TRUE, showCategory = 100,label_format = 30,font.size=8,color = 'pvalue',title = 'KEGG enrichment all')
dev.off()
# 绘制气泡图
pdf(file="0828.KO_enrichment_all_dot.pdf",width = 10,height = 7)
dotplot(result_all, showCategory = 30)


dev.off()


####部分富集结果####
KEGG_result <- data.frame(result_down)
KEGG_result <- KEGG_result[c(12:14,20:22,26:28),]
KEGG_result$GeneRatio_numeric <- sapply(KEGG_result$GeneRatio, convert_fraction_to_numeric)  
KEGG_result[order(KEGG_result$GeneRatio_numeric),]$Description
KEGG_result$Description<- factor(KEGG_result$Description,levels =KEGG_result[order(KEGG_result$GeneRatio_numeric),]$Description)


# 绘制条形图
pdf(file="0828.KO_enrichment_VB_bar.pdf",width = 10,height = 7)

ggplot(data = KEGG_result,
       aes(x = GeneRatio_numeric, y = Description, fill = p.adjust)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "KEGG") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()
