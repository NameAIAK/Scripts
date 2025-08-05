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
  id <- read_excel("HYXM_16S.id.xlsx",sheet = num)
  
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
  

  
  # 设置每个组的颜色
  group_colors <- c("AB" = "#BD3C29", "C" = "#78D3AC")
  # 首先设置比较的列表
  compare_list <- list(
    c("AB","C"))
  
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
  # taxonomy <- gsub("; s__", ";", taxonomy)
  rownames(count_tmp) <- taxonomy
  # 删除第一列OTU
  count_tmp <- count_tmp[,-1]
  # 删除最后一列taxonomy
  df <- count_tmp %>% select(-taxonomy)
  
  df_abun <- df
  
  df_abun <- df_abun[,c(meta_id$sample)]
  
  for (specie in c('Veillonella',
                   'Bacteroides',
                   'Bifidobacterium',
                   'Lachnospiraceae',
                   'Candidatus_Stoquefichus',
                   'Butyricimonas',
                   'Barnesiella',
                   'Megamonas'
                   )){
    print(specie)
    values1 <- df_abun[specie,]
    values1 <- t(values1)
    colnames(values1) <- c(specie)
    class(values1)
    
    metadata <- meta_id
    values1 <- as.data.frame(values1)
    values1$sample <- rownames(values1)
    
    plot_data <- merge(metadata, values1, by = "sample")
    
    
    # 使用ggplot2创建Shannon箱线图
    p_values1 <- ggplot(plot_data, aes(x = group, y =.data[[specie]], fill = group, colour = group)) +
      geom_boxplot(width = 0.5, alpha = 0.4, lwd = 1.15, outlier.shape = NA) +  # 调整箱的大小
      geom_jitter(width = 0.3, size = 3, alpha = 0.75) +  # 添加散点
      labs(y = sprintf("%s-%s",num,specie)) +
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
    ggsave(sprintf("%s-%s.pdf",num,specie), p_values1, height = 5, width = 5)
  }
  
  
}

####联合分析的样本信息####
for (num in c('1','2','3','4','5','6','7')){
  print(num)
  # num='1'
  id <- read_excel("HYXM_16S.id.xlsx",sheet = num)
  
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
  
  library(ggpubr)
  
  # 设置每个组的颜色
  # group_colors <- c("AB" = "#BD3C29", "C" = "#78D3AC")
  # 首先设置比较的列表
  compare_list <- list(
    c("AB","C"))
  
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
  
  df_abun <- df_abun[,c(meta_id$sample)]
  i=1
  fic<- list()
  for (specie in c('Veillonella',
                   'Bacteroides',
                   'Bifidobacterium',
                   'Candidatus_Stoquefichus',
                   'Butyricimonas',
                   'Barnesiella',
                   'Megamonas'
  )){
    print(specie)
    # specie='Bifidobacterium'
    values1 <- df_abun[specie,]
    values1 <- t(values1)
    colnames(values1) <- c(specie)
    class(values1)
    
    metadata <- meta_id
    values1 <- as.data.frame(values1)
    values1$sample <- rownames(values1)
    
    plot_data <- merge(metadata, values1, by = "sample")
    
    
    my_theme = theme_bw() +                                # 设置可选图表样式
      theme(plot.title = element_text(hjust = 0.5,size = 20),     # 设置标题文字居中和字体大小
            axis.text = element_text(size = 15,color = 'black'),  # 设置坐标轴标签文字的大小和颜色
            axis.title = element_text(size = 20))                 # 设置坐标轴标题的文字大小
    my_theme2 = my_theme + theme(legend.title = element_text(size = 20),
                                 legend.text = element_text(size = 12))
    
    fic[[i]]<- ggplot(plot_data, aes(x = group, y =.data[[specie]], fill = group, colour = group)) + 
      geom_boxplot() +
      guides(fill = 'none') + 
      labs(title = sprintf("%s-%s",num,specie)) + 
      my_theme2+ 
      scale_fill_manual(values = c('white','white')) +
      xlab(NULL)+
      ylab(NULL)
    i=i+1
    
  }
  ggsave(sprintf("%s.pdf",num),ggarrange(fic[[1]],fic[[2]],fic[[3]],fic[[4]],fic[[5]],fic[[6]],
            fic[[7]]),width = 13)
  
}
#####MY-CODE####
dat<- data.frame(1:43,1:43,1:43)
colnames(dat)<- c('sample','group','Level')
dat$sample<- colnames(dd)
dat$group<-c(rep('Treat',30),rep('Control',13))
dat$group<- factor(dat$group,levels = c('Control','Treat'))
fic<- list()
my_theme = theme_bw() +                                # 设置可选图表样式
  theme(plot.title = element_text(hjust = 0.5,size = 20),     # 设置标题文字居中和字体大小
        axis.text = element_text(size = 15,color = 'black'),  # 设置坐标轴标签文字的大小和颜色
        axis.title = element_text(size = 20))                 # 设置坐标轴标题的文字大小
my_theme2 = my_theme + theme(legend.title = element_text(size = 20),
                             legend.text = element_text(size = 12))

for (i in 1:dim(dd)[1]) {
  dat$Level<- as.numeric(dd[i,])
  fic[[i]]<- ggplot(dat, aes(group,Level,color=group,fill=group)) + 
    geom_boxplot() +
    guides(fill = 'none') + 
    labs(title = nam[i]) + 
    my_theme2+ 
    scale_fill_manual(values = c('white','white')) +
    
    xlab(NULL)
}

ggarrange(fic[[1]],fic[[2]],fic[[3]],fic[[4]],fic[[5]],fic[[6]],
          fic[[7]])

#####s-n######

for (specie in c('Veillonella',
                 'Bacteroides',
                 'Bifidobacterium',
                 'Candidatus_Stoquefichus',
                 'Butyricimonas',
                 'Barnesiella',
                 'Megamonas'
                 # 'annerostipes'
)){
  i=1
  fic<- list()
  for (num in c('1','2','3','4','5','6','7')){
   
    
    print(num)
    # num='1'
    id <- read_excel("HYXM_16S.id.xlsx",sheet = num)
    
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
    
    library(ggpubr)
    
    # 设置每个组的颜色
    # group_colors <- c("AB" = "#BD3C29", "C" = "#78D3AC")
    # 首先设置比较的列表
    compare_list <- list(
      c("Treat","Control"))
    
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
    
    df_abun <- df_abun[,c(meta_id$sample)]
    print(specie)
    # specie='Bifidobacterium'
    values1 <- df_abun[specie,]
    values1 <- t(values1)
    colnames(values1) <- c(specie)
    class(values1)
    
    metadata <- meta_id
    values1 <- as.data.frame(values1)
    values1$sample <- rownames(values1)
    
    plot_data <- merge(metadata, values1, by = "sample")
    
    
    my_theme = theme_bw() +                                # 设置可选图表样式
      theme(plot.title = element_text(hjust = 0.5,size = 20),     # 设置标题文字居中和字体大小
            axis.text = element_text(size = 15,color = 'black'),  # 设置坐标轴标签文字的大小和颜色
            axis.title = element_text(size = 20))                 # 设置坐标轴标题的文字大小
    my_theme2 = my_theme + theme(legend.title = element_text(size = 20),
                                 legend.text = element_text(size = 12))
    # 箱型图
    plot_data <- plot_data[plot_data[, 3] != 0, ]
    fic[[i]]<- ggplot(plot_data, aes(x = group, y =.data[[specie]], fill = group, colour = group)) +
      geom_boxplot() +
      # geom_point()+
      guides(fill = 'none') +
      labs(title = sprintf("%s-%s",num,specie)) +
      my_theme2+
      scale_fill_manual(values = c('white','white')) +
      xlab(NULL)+
      ylab(NULL)
    # # 散点图
    # plot_data$sample <- factor(plot_data$sample,levels = plot_data[order(plot_data$group),]$sample)
    # fic[[i]]<- ggplot(plot_data, aes(x = sample, y =.data[[specie]], fill = group, colour = group)) +
    #   # barplot() +
    #   geom_point()+
    #   guides(fill = 'none') +
    #   labs(title = sprintf("%s-%s",num,specie)) +
    #   my_theme2+
    #   scale_fill_manual(values = c('white','white')) +
    #   xlab(NULL)+
    #   ylab(NULL)+theme(axis.text.x=element_blank(),
    #                   axis.ticks.x=element_blank(),
    #                   )
    # # 柱状图
    # result <- plot_data %>%
    #   group_by(group) %>%
    #   summarise(
    #     mean_value = mean(.data[[specie]], na.rm = TRUE),
    #     sd_value = sd(.data[[specie]], na.rm = TRUE),
    # 
    #   )
    # fic[[i]]<- ggplot(result, aes(x=group, y=mean_value, , fill = group, colour = group)) +
    #     geom_bar(position=position_dodge(), stat="identity") +
    #     geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=.2)+
    #   # guides(fill = 'none') +
    #   labs(title = sprintf("%s-%s",num,specie)) + 
    #   my_theme2+ 
    #   # scale_fill_manual(values = c('white','white')) +
    #   xlab(NULL)+
    #   ylab(NULL)
    i=i+1
    
  }
  ggsave(sprintf("%s.pdf",specie),ggarrange(fic[[1]],fic[[2]],fic[[3]],fic[[4]],fic[[5]],fic[[6]],
                                         fic[[7]]),width = 13)

}

#####all#####
rm(list=ls())
i=1
fic<- list()
id <- read_excel("HYXM_16S.id.xlsx")
for (specie in c('Veillonella',
                 'Bacteroides',
                 'Bifidobacterium',
                 'Candidatus_Stoquefichus',
                 'Butyricimonas',
                 'Barnesiella',
                 'Megamonas'
                 # 'annerostipes'
)){
  # specie <- 'Megamonas'
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
  plot_data[plot_data[,4] == 0,4] <- 1e-9
  plot_data[,4] <- log(plot_data[,4],base = 10)
  plot_data[,4] <- scale(
    plot_data[,4]
  )
  fic[[i]]<- ggplot(plot_data, aes(x = group, y =.data[[specie]], fill = group, colour = group)) +
    geom_boxplot() +
    # geom_point()+
    guides(fill = 'none') +
    labs(title = sprintf("%s",specie)) +
    my_theme2+
    scale_fill_manual(values = c('white','white')) +
    xlab(NULL)+
    ylab(NULL)
  # # 散点图
  # plot_data$sample <- factor(plot_data$sample,levels = plot_data[order(plot_data$group),]$sample)
  # write.table(plot_data,file = sprintf("%s_abun.csv",specie), sep = ",", row.names = FALSE, quote = FALSE)
  # 
  # fic[[i]]<- ggplot(plot_data, aes(x = sample, y =.data[[specie]], fill = group, colour = group)) +
  #   # barplot() +
  #   geom_point()+
  #   guides(fill = 'none') +
  #   labs(title = sprintf("%s",specie)) +
  #   my_theme2+
  #   scale_fill_manual(values = c('white','white')) +
  #   # xlab(NULL)+
  #   ylab(NULL)
  # +theme(axis.text.x=element_blank(),
  #                   axis.ticks.x=element_blank(),
  #                   )
  # # 柱状图
  # result <- plot_data %>%
  #   group_by(group) %>%
  #   summarise(
  #     mean_value = mean(.data[[specie]], na.rm = TRUE),
  #     sd_value = sd(.data[[specie]], na.rm = TRUE),
  # 
  #   )
  # fic[[i]]<- ggplot(result, aes(x=group, y=mean_value, , fill = group, colour = group)) +
  #     geom_bar(position=position_dodge(), stat="identity") +
  #     # geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), width=.2)+
  #   # guides(fill = 'none') +
  #   labs(title = sprintf("%s",specie)) +
  #   my_theme2+
  #   # scale_fill_manual(values = c('white','white')) +
  #   xlab(NULL)+
  #   ylab(NULL)
  i=i+1
}
ggsave("allspecies.sample.point.pdf",ggarrange(fic[[1]],fic[[2]],fic[[3]],fic[[4]],fic[[5]],fic[[6]],
                                          fic[[7]]),width = 13)

