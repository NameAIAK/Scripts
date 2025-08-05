rm(list=ls())
library(tidyr)
library(dplyr)
library(readxl) 
library(data.table)
library(ggplot2)
library(ggpubr)
####联合分析的样本信息####
for (num in c('1','2','3','4','5','6','7')){
  print(num)
  # num='1'
  id <- read_excel("HYXM_16S.id.xlsx",sheet = num)
  # cor_abs <- 0.3
  
  
  PMB <- fread(sprintf("../PCancerMb/P-%s.csv",num))
  NMB <- fread(sprintf("../PCancerMb/N-%s.csv",num))
  
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
  PMB_NAME <- fread(sprintf('../Gred/r0.3/P-%s-差异离子鉴定.csv',num))
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
  # 
  rownames(mb_pos) <- rownames_mb_pos
  df_mb_pos <- t(mb_pos)
  
  
  NMB <- as.data.frame(NMB)
  colnames(NMB)[1] <- 'Compound'
  # NMB <- NMB[,c(1:50)]
  NMB_NAME <- fread(sprintf('../Gred/r0.3/N-%s-差异离子鉴定.csv',num))
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
  rownames_mb_neg <- mb_neg$name
  mb_neg <- mb_neg[,-1]
  
  rownames(mb_neg) <- rownames_mb_neg
  df_mb_neg <- t(mb_neg)
  
  Ergothioneine <- df_mb_neg[,c('Ergothioneine')]
  Ergothioneine <-as.data.frame(Ergothioneine)
  # colnames(Ergothioneine) <- 'Ergothioneine'
  # 获取前10个行名和数值
  # 获取数值最大的前10行名和对应的数值
  Ergothioneine <- Ergothioneine %>%
    arrange(desc(Ergothioneine)) 
  Ergothioneine$mb_id_N <- rownames(Ergothioneine)
  Ergothioneine <- merge(Ergothioneine,N_id,by='mb_id_N')
  write.csv(Ergothioneine,sprintf('%s.csv',num),row.names = TRUE)
  
}  


####abun_part####
rm(list=ls())
i=1
fic<- list()
id <- read_excel("./Ergothioneine_7_70.xlsx")
id2 <- read_excel("HYXM_16S.id.xlsx")

id_merge <-merge(id,id2,by='sample')
id <- id_merge[,c(6,2)]
colnames(id) <- c('sample','group')
table(id$group)
for (specie in c('Veillonella',
                 'Candidatus_Stoquefichus',
                 'Butyricimonas',
                 'Barnesiella'
                 # 'annerostipes'
)){
  # specie <- 'Veillonella'
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
  
  count_tmp <- data
  taxonomy <- sapply(strsplit(count_tmp$taxonomy, ";g__"), function(x) x[length(x)])
  # taxonomy <- sapply(strsplit(count_tmp$taxonomy, ";g__"), function(x) x[0])
  # taxonomy <- gsub("; s__", ";", taxonomy)
  rownames(count_tmp) <- taxonomy
  # 删除第一列OTU
  count_tmp <- count_tmp[,-1]
  # 删除最后一列taxonomy
  df <- count_tmp %>% select(-taxonomy)
  
  df_abun <- df
  
  # df_abun <- df_abun[,c(meta_id$sample)]
  print(specie)
  # specie='Bifidobacterium'
  values1 <- df_abun[specie,id$sample]
  values1 <- t(values1)
  # colnames(values1) <- c(specie)
  # class(values1)
  
  # metadata <- meta_id
  values1 <- as.data.frame(values1)
  sample <- rownames(values1)
  values1$sample <- sample
  
  plot_data <- merge(id, values1, by = "sample")
  
  
  my_theme = theme_bw() +                                # 设置可选图表样式
    theme(plot.title = element_text(hjust = 0.5,size = 20),     # 设置标题文字居中和字体大小
          axis.text = element_text(size = 15,color = 'black'),  # 设置坐标轴标签文字的大小和颜色
          axis.title = element_text(size = 20))                 # 设置坐标轴标题的文字大小
  my_theme2 = my_theme + theme(legend.title = element_text(size = 20),
                               legend.text = element_text(size = 12))
  # 箱型图
  # plot_data <- plot_data[plot_data[, 4] != 0, ]
  # plot_data[plot_data[,4] == 0,4] <- 1e-9
  # plot_data[,4] <- log(plot_data[,4],base = 10)
  # plot_data[,4] <- scale(
  #   plot_data[,4]
  # )
  # fic[[i]]<- ggplot(plot_data, aes(x = group, y =.data[[specie]], fill = group, colour = group)) +
  #   geom_boxplot() +
  #   # geom_point()+
  #   guides(fill = 'none') +
  #   labs(title = sprintf("%s",specie)) +
  #   my_theme2+
  #   scale_fill_manual(values = c('white','white')) +
  #   xlab(NULL)+
  #   ylab(NULL)
  # # 散点图
  # plot_data$sample <- factor(plot_data$sample,levels = plot_data[order(plot_data$group),]$sample)
  # # write.table(plot_data,file = sprintf("%s_abun.csv",specie), sep = ",", row.names = FALSE, quote = FALSE)
  # 
  # fic[[i]]<- ggplot(plot_data, aes(x = sample, y =.data[[specie]], fill = group, colour = group)) +
  #   # barplot() +
  #   geom_point()+
  #   guides(fill = 'none') +
  #   labs(title = sprintf("%s",specie)) +
  #   my_theme2+
  #   scale_fill_manual(values = c('white','white')) +
  #   # xlab(NULL)+
  #   ylab(NULL)+theme(axis.text.x=element_blank(),
  #                   axis.ticks.x=element_blank(),
  #                   )
  # 柱状图
  result <- plot_data %>%
    group_by(group) %>%
    summarise(
      mean_value = mean(.data[[specie]], na.rm = TRUE),
      sd_value = sd(.data[[specie]], na.rm = TRUE),

    )
  fic[[i]]<- ggplot(result, aes(x=group, y=mean_value, , fill = group, colour = group)) +
      geom_bar(position=position_dodge(), stat="identity") +
      # geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=.2)+
    # guides(fill = 'none') +
    labs(title = sprintf("%s",specie)) +
    my_theme2+
    # scale_fill_manual(values = c('white','white')) +
    xlab(NULL)+
    ylab(NULL)
  i=i+1
}
ggsave("E15.bar.pdf",ggarrange(fic[[1]],fic[[2]],fic[[3]],fic[[4]]),width = 13)


  


####abun_all####
rm(list=ls())
i=1
fic<- list()
# id <- read_excel("./Ergothioneine_7_70.xlsx")
id <- read_excel("HYXM_16S.id.xlsx",sheet = 'Sheet3')

# id_merge <-merge(id,id2,by='sample')
# id <- id_merge[,c(6,2)]
# colnames(id) <- c('sample','group')
# table(id$group)
for (specie in c('Veillonella',
                 'Candidatus_Stoquefichus',
                 'Butyricimonas',
                 'Barnesiella'
                 # 'annerostipes'
)){
  # specie <- 'Veillonella'
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
  
  count_tmp <- data
  taxonomy <- sapply(strsplit(count_tmp$taxonomy, ";g__"), function(x) x[length(x)])
  # taxonomy <- sapply(strsplit(count_tmp$taxonomy, ";g__"), function(x) x[0])
  # taxonomy <- gsub("; s__", ";", taxonomy)
  rownames(count_tmp) <- taxonomy
  # 删除第一列OTU
  count_tmp <- count_tmp[,-1]
  # 删除最后一列taxonomy
  df <- count_tmp %>% select(-taxonomy)
  
  df_abun <- df
  
  # df_abun <- df_abun[,c(meta_id$sample)]
  print(specie)
  # specie='Bifidobacterium'
  values1 <- df_abun[specie,id$sample]
  values1 <- t(values1)
  # colnames(values1) <- c(specie)
  # class(values1)
  
  # metadata <- meta_id
  values1 <- as.data.frame(values1)
  sample <- rownames(values1)
  values1$sample <- sample
  
  plot_data <- merge(id, values1, by = "sample")
  
  
  my_theme = theme_bw() +                                # 设置可选图表样式
    theme(plot.title = element_text(hjust = 0.5,size = 20),     # 设置标题文字居中和字体大小
          axis.text = element_text(size = 15,color = 'black'),  # 设置坐标轴标签文字的大小和颜色
          axis.title = element_text(size = 20))                 # 设置坐标轴标题的文字大小
  my_theme2 = my_theme + theme(legend.title = element_text(size = 20),
                               legend.text = element_text(size = 12))
  # 箱型图
  # plot_data <- plot_data[plot_data[, 4] != 0, ]
  # plot_data[plot_data[,4] == 0,4] <- 1e-9
  # plot_data[,4] <- log(plot_data[,4],base = 10)
  # plot_data[,4] <- scale(
  #   plot_data[,4]
  # )
  # fic[[i]]<- ggplot(plot_data, aes(x = group, y =.data[[specie]], fill = group, colour = group)) +
  #   geom_boxplot() +
  #   # geom_point()+
  #   guides(fill = 'none') +
  #   labs(title = sprintf("%s",specie)) +
  #   my_theme2+
  #   scale_fill_manual(values = c('white','white')) +
  #   xlab(NULL)+
  #   ylab(NULL)
  # 散点图
  plot_data$sample <- factor(plot_data$sample,levels = plot_data[order(plot_data$group),]$sample)
  # write.table(plot_data,file = sprintf("%s_abun.csv",specie), sep = ",", row.names = FALSE, quote = FALSE)

  fic[[i]]<- ggplot(plot_data, aes(x = sample, y =.data[[specie]], fill = group, colour = group)) +
    # barplot() +
    geom_point()+
    guides(fill = 'none') +
    labs(title = sprintf("%s",specie)) +
    my_theme2+
    scale_fill_manual(values = c('white','white')) +
    # xlab(NULL)+
    ylab(NULL)+theme(axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    )
  # # 柱状图
  # result <- plot_data %>%
  #   group_by(group) %>%
  #   summarise(
  #     mean_value = mean(.data[[specie]], na.rm = TRUE),
  #     sd_value = sd(.data[[specie]], na.rm = TRUE),
  #     
  #   )
  # fic[[i]]<- ggplot(result, aes(x=group, y=mean_value, , fill = group, colour = group)) +
  #   geom_bar(position=position_dodge(), stat="identity") +
  #   # geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=.2)+
  #   # guides(fill = 'none') +
  #   labs(title = sprintf("%s",specie)) +
  #   my_theme2+
  #   # scale_fill_manual(values = c('white','white')) +
  #   xlab(NULL)+
  #   ylab(NULL)
  i=i+1
}
ggsave("E15.point.all.pdf",ggarrange(fic[[1]],fic[[2]],fic[[3]],fic[[4]]),width = 13)

####abun_red####
rm(list=ls())
i=1
fic<- list()
# id <- read_excel("./Ergothioneine_7_70.xlsx")
id <- read_excel("HYXM_16S.id.xlsx",sheet = 'Sheet1')
id <- id[,c(1,2)]
# id_merge <-merge(id,id2,by='sample')
# id <- id_merge[,c(6,2)]
# colnames(id) <- c('sample','group')
# table(id$group)
for (specie in c('Veillonella',
                 'Candidatus_Stoquefichus',
                 'Butyricimonas',
                 'Barnesiella'
                 # 'annerostipes'
)){
  # specie <- 'Veillonella'
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
  
  count_tmp <- data
  taxonomy <- sapply(strsplit(count_tmp$taxonomy, ";g__"), function(x) x[length(x)])
  # taxonomy <- sapply(strsplit(count_tmp$taxonomy, ";g__"), function(x) x[0])
  # taxonomy <- gsub("; s__", ";", taxonomy)
  rownames(count_tmp) <- taxonomy
  # 删除第一列OTU
  count_tmp <- count_tmp[,-1]
  # 删除最后一列taxonomy
  df <- count_tmp %>% select(-taxonomy)
  
  df_abun <- df
  
  # df_abun <- df_abun[,c(meta_id$sample)]
  print(specie)
  # specie='Bifidobacterium'
  values1 <- df_abun[specie,id$sample]
  values1 <- t(values1)
  # colnames(values1) <- c(specie)
  # class(values1)
  
  # metadata <- meta_id
  values1 <- as.data.frame(values1)
  sample <- rownames(values1)
  values1$sample <- sample
  
  plot_data <- merge(id, values1, by = "sample")
  
  
  my_theme = theme_bw() +                                # 设置可选图表样式
    theme(plot.title = element_text(hjust = 0.5,size = 20),     # 设置标题文字居中和字体大小
          axis.text = element_text(size = 15,color = 'black'),  # 设置坐标轴标签文字的大小和颜色
          axis.title = element_text(size = 20))                 # 设置坐标轴标题的文字大小
  my_theme2 = my_theme + theme(legend.title = element_text(size = 20),
                               legend.text = element_text(size = 12))
  # 箱型图
  # plot_data <- plot_data[plot_data[, 4] != 0, ]
  # plot_data[plot_data[,4] == 0,4] <- 1e-9
  # plot_data[,4] <- log(plot_data[,4],base = 10)
  # plot_data[,4] <- scale(
  #   plot_data[,4]
  # )
  # fic[[i]]<- ggplot(plot_data, aes(x = group, y =.data[[specie]], fill = group, colour = group)) +
  #   geom_boxplot() +
  #   # geom_point()+
  #   guides(fill = 'none') +
  #   labs(title = sprintf("%s",specie)) +
  #   my_theme2+
  #   scale_fill_manual(values = c('white','white')) +
  #   xlab(NULL)+
  #   ylab(NULL)
  # 散点图
  plot_data$sample <- factor(plot_data$sample,levels = plot_data[order(plot_data$group),]$sample)
  # write.table(plot_data,file = sprintf("%s_abun.csv",specie), sep = ",", row.names = FALSE, quote = FALSE)
  
  fic[[i]]<- ggplot(plot_data, aes(x = sample, y =.data[[specie]], fill = group, colour = group)) +
    # barplot() +
    geom_point()+
    guides(fill = 'none') +
    labs(title = sprintf("%s",specie)) +
    my_theme2+
    scale_fill_manual(values = c('white','white')) +
    # xlab(NULL)+
    ylab(NULL)+theme(axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
    )
  # # 柱状图
  # result <- plot_data %>%
  #   group_by(group) %>%
  #   summarise(
  #     mean_value = mean(.data[[specie]], na.rm = TRUE),
  #     sd_value = sd(.data[[specie]], na.rm = TRUE),
  #     
  #   )
  # fic[[i]]<- ggplot(result, aes(x=group, y=mean_value, , fill = group, colour = group)) +
  #   geom_bar(position=position_dodge(), stat="identity") +
  #   # geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=.2)+
  #   # guides(fill = 'none') +
  #   labs(title = sprintf("%s",specie)) +
  #   my_theme2+
  #   # scale_fill_manual(values = c('white','white')) +
  #   xlab(NULL)+
  #   ylab(NULL)
  i=i+1
}
ggsave("E15.point.red.pdf",ggarrange(fic[[1]],fic[[2]],fic[[3]],fic[[4]]),width = 13)


id <- read_excel("./Ergothioneine_7_70.xlsx")
id2 <- read_excel("HYXM_16S.id.xlsx")

id_merge <-merge(id,id2,by='sample')
id <- id_merge[,c(6,2)]
colnames(id) <- c('sample','group')
table(id$group)
# for (specie in c('Veillonella',
#                  'Candidatus_Stoquefichus',
#                  'Butyricimonas',
#                  'Barnesiella'
#                  # 'annerostipes'
# )){
# specie <- 'Veillonella'
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

count_tmp <- data
taxonomy <- sapply(strsplit(count_tmp$taxonomy, ";g__"), function(x) x[length(x)])
# taxonomy <- sapply(strsplit(count_tmp$taxonomy, ";g__"), function(x) x[0])
# taxonomy <- gsub("; s__", ";", taxonomy)
rownames(count_tmp) <- taxonomy
# 删除第一列OTU
count_tmp <- count_tmp[,-1]
# 删除最后一列taxonomy
df <- count_tmp %>% select(-taxonomy)

df_abun <- df

# df_abun <- df_abun[,c(meta_id$sample)]
print(specie)
# specie='Bifidobacterium'
values1 <- df_abun[c('Veillonella','Candidatus_Stoquefichus',
                     'Butyricimonas','Barnesiella'),id$sample]
values1 <- t(values1)
# colnames(values1) <- c(specie)
# class(values1)

# metadata <- meta_id
values1 <- as.data.frame(values1)
sample <- rownames(values1)
values1$sample <- sample

plot_data <- merge(id, values1, by = "sample")
write.table(plot_data,file = "4species_abun.csv", sep = ",", row.names = FALSE, quote = FALSE)



####3mb####
for (num in c('1','2','3','4','5','6','7')){
  print(num)
  # num='1'
  id <- read_excel("HYXM_16S.id.xlsx",sheet = num)
  # cor_abs <- 0.3
  
  
  PMB <- fread(sprintf("../PCancerMb/P-%s.csv",num))
  NMB <- fread(sprintf("../PCancerMb/N-%s.csv",num))
  
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
  PMB_NAME <- fread(sprintf('../Gred/r0.3/P-%s-差异离子鉴定.csv',num))
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
  # 
  rownames(mb_pos) <- rownames_mb_pos
  df_mb_pos <- t(mb_pos)
  
  
  NMB <- as.data.frame(NMB)
  colnames(NMB)[1] <- 'Compound'
  # NMB <- NMB[,c(1:50)]
  NMB_NAME <- fread(sprintf('../Gred/r0.3/N-%s-差异离子鉴定.csv',num))
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
  rownames_mb_neg <- mb_neg$name
  mb_neg <- mb_neg[,-1]
  
  rownames(mb_neg) <- rownames_mb_neg
  df_mb_neg <- t(mb_neg)
  # 'Ergothioneine','Collettiside I','Cycloastragenol'
  df <- df_mb_neg[,c('Ergothioneine','Collettiside I','Cycloastragenol')]
  df <-as.data.frame(df)
  df <- t(df)
  # df$name <- rownames(df)
  write.csv(df,sprintf('%s-mb.csv',num),row.names = TRUE)
  
} 

####获取sample####
id1 <- read_excel("HYXM_16S.id.xlsx")
id2 <- read_excel("./E_species1-5/16+11sample_sabun.xlsx")
sample <- merge(id2,id1,by='meta_id')
sample <- sample[,c('meta_id','sample','group.y')]
colnames(sample) <- c('meta_id','sample','group')

mb3all <- read_excel('./3mb_abun/3mb-all.xlsx')
rownames_mb3all <- mb3all$...1
mb3all <- mb3all[,-1]
# rownames(mb3all) <- rownames_mb3all
mb3part <- mb3all[,sample$sample]
colnames(mb3part) <- sample$meta_id
rownames(mb3part) <- rownames_mb3all
class(mb3part)
mb3part_t <- t(mb3part)
meta_id <- rownames(mb3part_t)
mb3part_t <- as.data.frame(mb3part_t)
mb3part_t$meta_id <- meta_id

df <- merge(mb3part_t,sample,by='meta_id')

write.csv(df,sprintf('3mb_g.csv'),row.names = TRUE)
