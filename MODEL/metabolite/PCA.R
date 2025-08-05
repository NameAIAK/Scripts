#将R自带的范例数据集iris储存为变量data;
data<-iris
head(data)
#对原数据进行z-score归一化；
dt<-as.matrix(scale(data[,1:4]))
head(dt)

#计算相关系数矩阵；
rm1<-cor(dt)
rm1
# 求解特征值和相应的特征向量
rs1<-eigen(rm1)
rs1
#提取结果中的特征值，即各主成分的方差；
val <- rs1$values
#换算成标准差(Standard deviation);
(Standard_deviation <- sqrt(val))
#计算方差贡献率和累积贡献率；
(Proportion_of_Variance <- val/sum(val))
(Cumulative_Proportion <- cumsum(Proportion_of_Variance))
#碎石图绘制;
par(mar=c(6,6,2,2))
plot(rs1$values,type="b",
cex=2,
cex.lab=2,
cex.axis=2,
lty=2,
lwd=2,
xlab = "主成分编号",
ylab="特征值(主成分方差)")

#提取结果中的特征向量(也称为Loadings,载荷矩阵)；
(U<-as.matrix(rs1$vectors))
#进行矩阵乘法，获得PC score；
PC <-dt %*% U
colnames(PC) <- c("PC1","PC2","PC3","PC4")
head(PC)

#将iris数据集的第5列数据合并进来；
df<-data.frame(PC,iris$Species)
head(df)

#载入ggplot2包；
library(ggplot2)
#提取主成分的方差贡献率，生成坐标轴标题；
xlab<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
ylab<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")
#绘制散点图并添加置信椭圆；
p1<-ggplot(data = df,aes(x=PC1,y=PC2,color=iris.Species))+
stat_ellipse(aes(fill=iris.Species),
type ="norm", geom ="polygon",alpha=0.2,color=NA)+
geom_point()+labs(x=xlab,y=ylab,color="")+
guides(fill=F)
p1

#载入scatterplot3d包；
library(scatterplot3d)
color = c(rep('purple',50),rep('orange',50),rep('blue',50))
scatterplot3d(df2[,1:3],color=color,
pch = 16,angle=30,
box=T,type="p",
lty.hide=2,lty.grid = 2)
legend("topleft",c('Setosa','Versicolor','Virginica'),
fill=c('purple','orange','blue'),box.col=NA)

# 使用现成函数完成主成分分析
# 1.prcomp()函数
# scale. = TRUE表示分析前对数据进行归一化；
com1 <- prcomp(data[,1:4], center = TRUE,scale. = TRUE)
summary(com1)

# 2.princomp()函数
com2 <- princomp(dt,cor = T)
summary(com2)
com3 <- princomp(dt)
summary(com3)

# 3.PCA结果可视化
#提取PC score；
df1<-com1$x
head(df1)
#将iris数据集的第5列数据合并进来；
df1<-data.frame(df1,iris$Species)
head(df1)

#提取主成分的方差贡献率，生成坐标轴标题；
summ<-summary(com1)
xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
p2<-ggplot(data = df1,aes(x=PC1,y=PC2,color=iris.Species))+
stat_ellipse(aes(fill=iris.Species),
type = "norm", geom ="polygon",alpha=0.2,color=NA)+
geom_point()+labs(x=xlab,y=ylab,color="")+
guides(fill=F)
p2+scale_fill_manual(values = c("purple","orange","blue"))+
scale_colour_manual(values = c("purple","orange","blue"))

