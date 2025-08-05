setwd('D:\\data\\李丛舒\\')
library(org.Mm.eg.db) ##加载人类
keytypes(org.Mm.eg.db)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(clusterProfiler)

d1<- read.table('dat\\C_20_Clean_count.txt')
dat<- list()
all_files<-list.files('./dat')
for (i in 1:length(all_files)) {
  temp<- read.table(paste0('dat\\',all_files[i]))
  colnames(temp)<- temp[1,]
  temp<- temp[-1,]
  dat[[i]]<- temp
  print(i)
}
data<- dat[[1]]
for (i in 2:length(dat)) {
  data<- cbind(data,dat[[i]])
  
}
colnames(data)
rownames(data)<- data[,1]
data<- data[,-seq(1,24,2)]
colnames(data)

colnames(data)<- c('C20','C37','C90','C91',
                   'M21','M34','M52','M68',
                   'RH14','RH27','RH66','RH98')
#data<- data[,c('B12','B15','B362','B36','B4','C19','C27','C32','C34','C8','M2','M29','M37','M5','M7')]
gene<- rownames(data)
for (i in 1:length(gene)) {
  temp<- strsplit(gene[i],'[.]')[[1]][1]
  gene[i]<- temp
}
rownames(data)<- gene
gs<- bitr(rownames(data), fromType = "ENSEMBL", toType="SYMBOL",OrgDb = org.Mm.eg.db)
table(duplicated(gs$ENSEMBL))
table(duplicated(gs$SYMBOL))
gs<- gs[!duplicated(gs$ENSEMBL),]
gs<- gs[!duplicated(gs$SYMBOL),]
data<- data[gs$ENSEMBL,]
rownames(data)<-gs$SYMBOL
###TMM标化
library(edgeR)
head(data)
for (i in 1:dim(data)[2]) {
  data[,i]<- as.numeric(data[,i])
  
}
group<- factor(c(rep('C',3),rep('M',3),rep('RH',3)))
#group<- factor(c(rep('B',5),rep('C',5),rep('M',5)))
y<- DGEList(counts = data,group = group)
y<- calcNormFactors(y, method = "TMM")
y$samples
norm_counts<- cpm(y, normalized.lib.sizes = TRUE)
head(norm_counts)
write.csv(data,file = 'count_2.csv')
write.csv(norm_counts,file = 'TMM_2.csv')


#差异表达基因分析
#首先根据分组信息构建试验设计矩阵，分组信息中一定要是对照组在前，处理组在后
colnames(data)
data<- data[,-c(4,8,11)]
dd<- data[,c(1:3,7:9)]
colnames(dd)
group2<-  rep(c('C', 'RH'), each = 3)
dgelist <- DGEList(counts = dd, group = group2)

#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]

#（3）标准化，以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#差异表达基因分析

design <- model.matrix(~group2)

#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)

#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))

write.table(lrt, 'RHvsC_2.txt', sep = '\t', col.names = NA, quote = FALSE)

#####
df2<- data.frame(1:20,1:20,1:20,1:20,1:20,1:20,1:20,1:20)
colnames(df2)<- c('names','group','Tnf','Il6','Cxcl1','Il1b','Cxcl2','Nfe2l2')
df2[,1]<- colnames(norm_counts)
df2[,2]<- c(rep('A',5),rep('B',5),rep('C',5),rep('M',5))
df2[,3]<- as.numeric(norm_counts['Tnf',])
df2[,4]<- as.numeric(norm_counts['Il6',])
df2[,5]<- as.numeric(norm_counts['Cxcl1',])
df2[,6]<- as.numeric(norm_counts['Il1b',])
df2[,7]<- as.numeric(norm_counts['Cxcl2',])
df2[,8]<- as.numeric(norm_counts['Nfe2l2',])

table(df2$group)
df2$group<- factor(df2$group,levels = c('C','M','A','B'))
p6 <- ggplot(df2,aes(x=group,y=Nfe2l2))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
  geom_boxplot(aes(fill=group), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(),#坐标轴的线设为显示
        plot.title = element_text(size=14))+#图例位置
  # scale_fill_manual(values=c("#ffc000","#a68dc8","blue"))+#指定颜色
  geom_jitter(width = 0.2)+#添加抖动点
  geom_signif(comparisons = list(c("A","M"),
                                 c("B","M"),
                                 c("C","M")),# 设置需要比较的组
              map_signif_level = F, #是否使用星号显示
              test = t.test, ##计算方法
            
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0)),#横线下方的竖线设置
              size=0.8,color="black")+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45)+ # 可选值有 0，45，90，270
  scale_fill_prism(palette = "candy_bright")+
  theme(legend.position = 'none')#去除图例
p6

ggarrange(p1,p2,p3,p4,p5,p6)

y_position = c(3,3.5,3.25)#图中横线位置 设置


df2<- df2[-17,]
df2<- df2[-c(1,7,15),]
p6<- ggplot(df2,aes(x=group,y=Nfe2l2))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
  geom_boxplot(aes(fill=group), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(),#坐标轴的线设为显示
        plot.title = element_text(size=14))+#图例位置
  # scale_fill_manual(values=c("#ffc000","#a68dc8","blue"))+#指定颜色
  geom_jitter(width = 0.2)+#添加抖动点
  geom_signif(comparisons = list(c("A","M"),
                                 c("B","M"),
                                 c("C","M")),# 设置需要比较的组
              map_signif_level = F, #是否使用星号显示
              test = t.test, ##计算方法
              y_position = c(28,27.5,27),
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0)),#横线下方的竖线设置
              size=0.8,color="black")+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45)+ # 可选值有 0，45，90，270
  scale_fill_prism(palette = "candy_bright")+
  theme(legend.position = 'none')#去除图例
p6










