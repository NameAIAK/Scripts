norm_counts<- norm_counts[,-c(4,8,11)]

PCAmat_t<-t(norm_counts)
PCAmat_t[1:5,1:3]
library(FactoMineR)
library(ggplot2)
library(ggrepel)
#样本中基因表达值的 PCA 分析
gene.pca <- PCA(PCAmat_t, ncp = 3, scale.unit = TRUE, graph = FALSE)

#提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:3])
pca_sample$Sample=row.names(pca_sample)
#提取 PCA 前两轴的贡献度(22.46,19.13,12.62)
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )

head(pca_sample)
rownames(pca_sample)
pca_sample$group<- c(c(rep('C',3),rep('M',3),rep('RH',3)))
#pca_sample$group<- c(c(rep('B',5),rep('C',5),rep('M',5)))

library(ggplot2)
p <- ggplot(data = pca_sample, aes(x =Dim.1 , y = Dim.2)) +
  geom_point(aes(color = group), size = 4.5) +  #根据样本坐标绘制二维散点图
   #自定义颜色
  scale_color_manual(values = c("#1E90FF","#E41A1C","green"))+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  #将 PCA 轴贡献度添加到坐标轴标题中

p

library(ggrepel)
p <- p + 
  geom_text_repel(aes(label = Sample), size = 2.5, show.legend = FALSE, 
                  box.padding = unit(0.2, 'lines'))

p

p + stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)

p + stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE) 

####3D
library(scatterplot3d)
dfGROUP<- data.frame(samples=colnames(norm_counts),group=colnames(norm_counts))
dfGROUP$group<-  c(c(rep('C',3),rep('M',3),rep('RH',3)))
#dfGROUP$group<-  c(c(rep('B',5),rep('C',5),rep('M',5)))
# PCA计算
pca_result <- prcomp(t(norm_counts),
                     scale=F  # 一个逻辑值，指示在进行分析之前是否应该将变量缩放到具有单位方差
)
pca_result$x<-data.frame(pca_result$x)


colors <-c(rep("#1E90FF",3),rep("#E41A1C",3),rep('green',3))
#colors <-c(rep("#1E90FF",5),rep("#E41A1C",5),rep('green',5))


# 绘图
s3d <- scatterplot3d(pca_result$x[,1:3],
                     pch = 16,       # 点形状
                      color = colors,  # 点颜色
                     cex.symbols = 2 # 点大小
)
# 设置图例
legend("top",
       legend = c('C','M','RH'),
       col =c("#1E90FF","#E41A1C","green") ,
       pch = 16,
       inset = -0.1,
       xpd = TRUE,
       horiz = TRUE)
#legend("top",legend = c('B','C','M'),col =c("#1E90FF","#E41A1C","green") ,pch = 16,inset = -0.1,xpd = TRUE,horiz = TRUE)
# 设置文字标注
text(s3d$xyz.convert(pca_result$x[,c(1,2,3)]),
     labels = row.names(pca_result$x),
     cex = 0.6,col = "black")


####3D_2
# 导入包
library(rgl)
options(stringsAsFactors = FALSE) #禁止chr转成factor
PCAmat_t[1:3,1:3]

## pca分析
pca <- prcomp(t(PCAmat_t))

# 01.导入分类数据
dfGROUP<- data.frame(samples=rownames(PCAmat_t),group=rownames(PCAmat_t))
dfGROUP$group<- c(c(rep('C',3),rep('M',3),rep('RH',3)))
#dfGROUP$group<- c(c(rep('B',5),rep('C',5),rep('M',5)))


## 
# 准备颜色
colors <-c(rep("#1E90FF",3),rep("#E41A1C",3),rep('green',3))
#colors <-c(rep("#1E90FF",5),rep("#E41A1C",5),rep('green',5))
#ctrl,LPS,YJSTW16h,YJCTW16h,PJHASTW16h,PJHACTW16h,EZSTW16h,EZCTW16h

plot3d(pca_result$x[,1:3], # 取前三个主成分
       xlab="PC1", ylab="PC2", zlab="PC3", 
       col=colors, # 按groups填充颜色
       type="s", # 画球，'p' for points, 's' for spheres, 'l' for lines, 'h' for line segments 
       size=2, #球的大小
       lwd=2, box=T)

#####################3d_3
library(pca3d)
data(metabo)

rotation<- pca_result$rotation[,1:3]
colnames(rotation)
PC1<- rotation[,'PC1']
names(PC1)<- rownames(rotation)
PC1[1:5]
sort(PC1,decreasing = T)[1:10]

##
PC2<- rotation[,'PC2']
names(PC2)<- rownames(rotation)
PC2[1:5]
sort(PC2,decreasing = T)[1:10]

##
PC3<- rotation[,'PC3']
names(PC3)<- rownames(rotation)
PC3[1:5]
sort(PC3,decreasing = T)[1:10]
