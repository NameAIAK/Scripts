  # ## 解包安装
  # # 下载安装包
  # wget -c ftp://download.nmdc.cn/tools//conda/lefse.tar.gz
  # 
  # # 解压到你的软件环境
  # tar -xvzf lefse.tar.gz -C xxx/envs/lefse
  # 
  # # conda启动环境
  # conda activate lefse
  # 
  # ## conda一键安装（经常出问题）
  # conda create -n lefse lefse -c bioconda -y
  # 
  # 
  # # 格式转换为lefse内部格式
  # lefse-format_input.py metaphlan4/lefse.txt \
  #   temp/input.in -c 2 -u 1 -o 1000000
  # 
  # 
  #   # 运行lefse
  # run_lefse.py input.in input.res
  # 
  # # 绘制物种分类层级树
  # lefse-plot_cladogram.py input.res \
  # lefse_cladogram.pdf --format pdf
  # 
  # # 绘制所有差异物种LDA得分柱状图
  # lefse-plot_res.py input.res \
  # lefse_res.pdf --format pdf
  # 
  # # 查看显著差异物种，按丰度排序
  # grep -v '-' input.res | sort -k3,3n 
  # 
  # # 批量绘制所有差异物种LDA得分柱状图
  # lefse-plot_features.py -f diff \
  #   --archive none --format pdf \
  #   input.in input.res \
  #   lefse
  # 加载tidyr包  
  # install.packages('tidyr')
####开始####
rm(list=ls())
library(tidyr)
library(dplyr)

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
  
otu <- data
original_colnames <- colnames(otu) 
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
colnames(otu) <- cleaned_colnames 
otu <- otu %>% select(-taxonomy)
rownames(otu) <- otu$`#OTU ID`
otu <- otu[,-1]


# 使用separate()函数分割combined_data列，并生成6列新数据  
df_split <- separate(data, col = taxonomy, into = paste0("col", 1:7), sep = ";")  
rownames(df_split) <- df_split$`#OTU ID`
df_split <- df_split[,-1]
last_7_cols <- df_split[, -(1:(ncol(df_split) - 7))] 
colnames(last_7_cols) <- c('Kingdom','Phylum',
                        'Class','Order',
                        'Family','Genus',
                        'Species')
tax <- last_7_cols

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
map$sample <- rownames(map)
colnames(map) <- c('group','sample')
# sample <- map



##加载R包

library(microeco) # Microbial Community Ecology Data Analysis
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(magrittr) # A Forward-Pipe Operator for R
library(ggtree)

# install.packages('ggtree')
#统一分类信息
tax %<>% tidy_taxonomy

## 构造microtable
df <- microtable$new(sample_table = map,
                     otu_table = otu,
                     tax_table = tax,
                     auto_tidy = F)
# df

# ####数据处理步骤，根据个人数据进行操作
# ##去除不属于非古菌和细菌的OTU
# df$tax_table %<>% subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
# df$tax_table %<>% .[grepl("Bacteria|Archaea", .$Kingdom), ]
# df
# 
# ##去除“线粒体”和“叶绿体”污染
df$filter_pollution(taxa = c("mitochondria", "chloroplast"))
# df
# 
# ##统一各数据的样本和OTU信息
df$tidy_dataset()
# df
# 
# ##检查序列号
df$sample_sums() %>% range
# 
# ##重采样以减少测序深度对多样性测量的影响，使每个样本的序列号相等。
df$rarefy_samples(sample.size = 10000)
df$sample_sums() %>% range
# 3、lefse分析：
# 主要基于microeco包中的trans_diff$new()函数实现，设置其中的方法为"lefse"并设置好分组即可，运行后命令行就会提示具体结果：
df$sample_table$group %<>% factor(., levels = c("CONTROL","OSAHS"))
str(df$sample_table)
# 
df$tax_table %<>% tidy_taxonomy
df <- na.omit(df) 
# df$cal_abund()

lefse <- trans_diff$new(dataset = df, #数据
                        method = "lefse", #方法
                        group = "group",#分组
                        lefse_subgroup = NULL,
                        p_adjust_method = "none"
)
# trans_diff$plot_diff_bar(
# color_values = RColorBrewer::brewer.pal(8, "Dark2"),
# color_group_map = FALSE,
# use_number = 1:10,
# threshold = NULL,
# select_group = NULL,
# keep_full_name = FALSE,
# keep_prefix = TRUE,
# group_order = NULL,
# group_aggre = TRUE,
# group_two_sep = TRUE,
# coord_flip = TRUE,
# add_sig = FALSE,
# add_sig_increase = 0.1,
# add_sig_text_size = 5,
# xtext_angle = 45,
# xtext_size = 10,
# axis_text_y = 12,
# heatmap_cell = "P.unadj",
# heatmap_sig = "Significance",
# heatmap_x = "Factors",
# heatmap_y = "Taxa",
# heatmap_lab_fill = "P value",
# ...
# )

# 取LDA>=2
p_LDA <- lefse$plot_diff_bar(threshold = 2.5,#设定LDA SCORE显示阈值
                    width = 0.8,#柱子宽度
                    group_order = c( "OSAHS", "CONTROL"),
                    color_values = RColorBrewer::brewer.pal(8, "Set1"))#分组顺序

p_LDA
ggsave(
  plot = p_LDA,
  "0828.p_LDA.pdf",
  height = 8,
  width = 6
)

# trans_diff$plot_diff_abund(
#   use_number = 1:20,
#   color_values = RColorBrewer::brewer.pal(8, "Dark2"),
#   select_group = NULL,
#   select_taxa = NULL,
#   simplify_names = TRUE,
#   keep_prefix = TRUE,
#   group_order = NULL,
#   barwidth = 0.9,
#   use_se = TRUE,
#   add_sig = FALSE,
#   add_sig_label = "Significance",
#   add_sig_label_color = "black",
#   add_sig_tip_length = 0.01,
#   y_start = 1.01,
#   y_increase = 0.05,
#   text_y_size = 10,
#   coord_flip = TRUE,
#   xtext_angle = 45,
#   ...
# )
# 默认取LDA前20
p_abund <- lefse$plot_diff_abund(group_order = c("OSAHS", "CONTROL"),
                                 use_number = 1:20,
                      add_sig = T,
                      color_values = RColorBrewer::brewer.pal(8, "Set1"))#是否显示差异
p_abund
ggsave(
  plot = p_abund,
  "0828.p_abund1-20.pdf",
  height = 6,
  width = 8
)

# 
# p_cladogram <- lefse$plot_diff_cladogram(use_taxa_num = 200,#显示丰度最高的40个分类群
#                           use_feature_num = 50, #显示30个差异特征
#                           clade_label_level = 5, #用字母标记标签的分类层级，5表示目，4表示科。
#                           # group_order = c("O", "C"),#组间排序
#                           color = RColorBrewer::brewer.pal(8, "Dark2"),#颜色
#                           filter_taxa = NULL,#过滤掉相对丰度低于0.0001的分类单元。
#                           select_show_labels = NULL,
#                           only_select_show = TRUE,# 设置为TRUE，可以只展示select_show_labels选择的分类单元标签
#                           sep = "|",#识别的辨识字符间隔
#                           branch_size = 0.5,
#                           alpha = 0.2,
#                           clade_label_size = 1.5)

# trans_diff$plot_diff_cladogram(
#   color = RColorBrewer::brewer.pal(8, "Dark2"),
#   group_order = NULL,
#   use_taxa_num = 200,
#   filter_taxa = NULL,
#   use_feature_num = NULL,
#   clade_label_level = 4,
#   select_show_labels = NULL,
#   only_select_show = FALSE,
#   sep = "|",
#   branch_size = 0.2,
#   alpha = 0.2,
#   clade_label_size = 2,
#   clade_label_size_add = 5,
#   clade_label_size_log = exp(1),
#   node_size_scale = 1,
#   node_size_offset = 1,
#   annotation_shape = 22,
#   annotation_shape_size = 5
# )
library(ggtree)
p_cladogram <-lefse$plot_diff_cladogram( 
                     clade_label_level = 5,
                     color = RColorBrewer::brewer.pal(8, "Set1"),
                     group_order = c("OSAHS", "CONTROL"))

p_cladogram
ggsave(
  plot = p_cladogram,
  "0828.p_cladogram.pdf",
  height = 10,
  width = 14
)



####随机森林分析####
# install.packages('randomForest')
random_forest_all <- trans_diff$new(dataset = df, method = "rf", 
                                group = "group", 
                                taxa_level = "all",
                                p_adjust_method = "none")#指定分类等级，如果使用所有分类群，可更改为"all"
# 绘制MeanDecreaseGini条形图
p_rfgini_all <- random_forest_all$plot_diff_bar(use_number = 1:20, 
                                  group_order = c("OSAHS", "CONTROL"),
                                  color_values = RColorBrewer::brewer.pal(8, "Set1"))
ggsave(
  plot = p_rfgini_all,
  "0828.p_rfgini.all.pdf",
  height = 6,
  width = 6
)
# 使用p1中相同的分类群绘制丰度图

p_rfabun_all <- random_forest_all$plot_diff_abund(use_number = 1:20,group_order = c("OSAHS", "CONTROL"),
                                    add_sig = T,
                                    color_values = RColorBrewer::brewer.pal(8, "Set1"))
ggsave(
  plot = p_rfabun_all,
  "0828.p_rfabun1-20.all.pdf",
  height = 6,
  width = 8
)
# 拼图
# p1 <- p1 + theme(legend.position = "none")
# p2 <- p2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
# gridExtra::grid.arrange(p1, p2, ncol = 2, nrow = 1, widths = c(1, 1.2))

####rf-genus####
random_forest <- trans_diff$new(dataset = df, method = "rf", 
                                group = "group", 
                                taxa_level = "Genus",
                                p_adjust_method = "none")#指定分类等级，如果使用所有分类群，可更改为"all"
# 绘制MeanDecreaseGini条形图
p_rfgini <- random_forest$plot_diff_bar(use_number = 1:20, 
                                        group_order = c("OSAHS", "CONTROL"),
                                        color_values = RColorBrewer::brewer.pal(8, "Set1"))
ggsave(
  plot = p_rfgini,
  "0828.p_rfgini.genus.pdf",
  height = 6,
  width = 6
)
# 使用p1中相同的分类群绘制丰度图

p_rfabun <- random_forest$plot_diff_abund(use_number = 1:20,group_order = c("OSAHS", "CONTROL"),
                                          add_sig = T,
                                          color_values = RColorBrewer::brewer.pal(8, "Set1"))
ggsave(
  plot = p_rfabun,
  "0828.p_rfabun1-20.genus.pdf",
  height = 6,
  width = 8
)
  
  
