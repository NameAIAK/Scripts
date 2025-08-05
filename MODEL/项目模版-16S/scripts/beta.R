# install.packages("ggplot2")  
#设置工作目录
# setwd("./")
rm(list=ls())
###########获取beta_diversity############
# install.packages('ggprism')
# 展示各组样本的微生物组成差异及组间差异
library(vegan)
library(ggplot2)
library(dplyr)
library(ggprism)
# install.packages("ggExtra")
library(ggExtra)
otu_table  <-  read.delim("./HYXM_16S.G.count.tsv.all.tmp", header=T, sep="\t", row.names=1, stringsAsFactors = FALSE)
otu_table <- otu_table %>% select(-taxonomy)
original_colnames <- colnames(otu_table) 
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
colnames(otu_table) <- cleaned_colnames  

# 使用vegan包中的Anosim函数、MRPP函数及Adonis函数进行进一步组间差异分析
#数据处理及PCoA分析
otu <- t(otu_table)
group <- read.table("HYXM_16S.red.group.txt", sep='\t')
colnames(group) <- c("samples","group")
otu <- otu[group$samples,]

otu.distance <- vegdist(otu)
PCoA <- cmdscale (otu.distance,eig=TRUE)
pc12 <- PCoA$points[,1:2]
pc <- round(PCoA$eig/sum(PCoA$eig)*100,digits=2)#解释度
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)

df <- merge(pc12,group,by="samples")
head(df)
colnames(df) <- c( 'samples','PCoA1','PCoA2','group')


# 进行置换多元（因素）方差分析（PERMANOVA）
nrow(otu.distance)  # 检查otu.distance的行数  
nrow(group)         # 检查group的行数
summary(otu.distance)  # 查看是否有NA  
summary(group)         # 查看是否有NA 
# 基于bray-curtis距离进行
dune.div <- adonis2(otu.distance ~ group, data = df, permutations = 999, method="bray")
dune.div
dune_adonis <- paste0("adonis R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)
dune_adonis

######pcoa######
p_pcoa <- ggplot(df,aes(x=PCoA1, y=PCoA2,color=group,shape=group))+#指定数据、X轴、Y轴
  geom_point(size=3)+labs(x=paste("PCo 1 (", pc[1], "%)", sep=""),
                          y=paste("PCo 2 (", pc[2], "%)", sep=""),
                          caption = dune_adonis)  +
  theme(legend.position = c(0.9,0.1),
        legend.title = element_blank(),
        panel.grid = element_blank(),plot.title=element_text(hjust=0.5),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        axis.text = element_text(color = "black",size=10))

p_pcoa <- p_pcoa+ stat_ellipse(data=df, 
                     geom = "polygon", 
                     level=0.9, 
                     linetype = 2, 
                     linewidth=0.5, 
                     aes(fill=group), 
                     alpha=0.3, 
                     show.legend = T) +
  scale_fill_manual(values = c("#f2ccac","#a1d5b9","#e1abbc")) 
p_pcoa

pdf(file="RHYXM_16S.pcoa.beta.pdf", width=5, height=5) 
## pdf格式
# pdf('marginal.pdf',width = 5,height = 5)
ggMarginal(
  p_pcoa,
  type =c('density'),
  margins = 'both',
  size = 3.5,
  groupColour = F,
  groupFill = T
) 
dev.off()
# Adonis,P>0.05,表示不同组的样品之间不存在显著差异

#设置随机种子
# set.seed(123)
# #基于bray-curtis距离进行PERMANOVA分析
# adonis <-  adonis2(dune ~ Management, data = dune.env, permutations = 999, method = "bray")
# #基于bray-curtis距离进行anosim分析
# anosim = anosim(dune, dune.env$Management, permutations = 999, distance = "bray")


#使用vegan包中的anosim函数进行anosim分析
# 当R趋向于1时，说明组间差异大于组内差异
# 当R=0时，说明组间没有差异，即分组无效，不同分组之间没有差异。
# 当R趋向于-1时，说明组间差异小于组内差异。
df_anosim <- anosim(otu.distance,df$group,permutations = 999)#数据也可以是原始otu数据
#df_anosim <- anosim(otu,df$group,permutations = 999)
#整理出作图数据
df1<-data.frame(
  x=df_anosim$class.vec,
  y=df_anosim$dis.rank
)
#绘图
p_anosim <- ggplot(df1,aes(x=x,y=y))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条先不会被箱体覆盖
  geom_boxplot(aes(fill=x), 
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        legend.position="none",plot.title = element_text(size=14))+
  scale_fill_manual(values=c("#1597A5","#FFC24B","#FEB3AE","red"))+ #指定颜色
  ggtitle("Bray-Curtis Anosim")+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain",
              base_family = "serif", 
              base_size = 14,  
              base_line_size = 0.8, 
              axis_text_angle = 45)+
  theme(legend.position = 'none')+
  labs(x = paste("R=",df_anosim$statistic,", ","p=", df_anosim$signif),
       y = "Rank of Distance (Bray_Curtis)")
p_anosim
ggsave("RHYXM_16S.anosim.beta.pdf", p_anosim, height = 5, width = 5)

# # 通过查看其中的A值与Pvalue值即可
# MRPP <- mrpp(otu.distance,df$group,permutations = 999)
# MRPP
# 
# 
# # 通过查看其中的R^2值及P值检验分组是否合理，其中R^2在0-1之间，越大则代表分组指标对差异的稀释度越高。
Adonis <- adonis2(otu.distance~group,data=df,
                  distance = "bray",
                  permutations = 999)
Adonis
# 
# 
# 
# 
# hc <- hclust(as.dist(otu.distance))  
# plot(hc)

####NMDS####
# NMDS排序，定义3个维度
nmds <- metaMDS(otu.distance, k = 3)

# 应力函数值，一般不大于 0.2 为合理
stress=round(nmds$stress, 3)
stress

# 样方得分
nmds_dis<- data.frame(nmds$points)
# write.table(nmds_dis_site,'nmds_dis_site.txt', sep = '\t', col.names = NA, quote = FALSE)
#####NMDS 评估######
#NMDS 评估，拟合 R2 越大越合理；
#检查观测值非相似性与排序距离之间的关系——没有点分布在线段较远位置表示该数据可以使用NMDS分析
stressplot(nmds)

# 提取列名，便于后面操作。
nmds_dis$Sample_ID <- rownames(nmds_dis)

# 和group合并
nmds_result <- cbind(nmds_dis,group)
p_v <- round(Adonis$`Pr(>F)`, 2)[1]
p_sig <- ifelse(p_v < 0.001, "***", ifelse(p_v < 0.01, "**", ifelse(p_v < 0.05, "*", "")))
adonis_text <- paste(paste("Adonis_p  =", round(Adonis$`Pr(>F)`, 2)), p_sig)[1]

p1 <- ggplot(data = nmds_result, aes(x=MDS1, y=MDS2, color=group)) +
  geom_point(aes(color=group,shape=group),size=5)+
  labs(x=paste("NMDS 1"),
       y=paste("NMDS 2"),
       caption=paste(paste('Stress =', stress),adonis_text))  +# 也可用title、caption
  scale_colour_manual(values = c("#314826", "#ca7442", "#5a558b"))+
  theme(legend.position = c(0.2,0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),plot.title=element_text(hjust=0),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        axis.text = element_text(color = "black",size=10))+
  geom_hline(aes(yintercept=0), colour="#BEBEBE", linetype="dashed")+
  geom_vline(aes(xintercept=0), colour="#BEBEBE", linetype="dashed")

p1
# ggsave('NMDS1.png',p1, width = 5, height = 5)


# 添加置信椭圆
p2 = p1 + stat_ellipse(data=nmds_result, 
                       geom = "polygon", 
                       level=0.9, 
                       linetype = 2, 
                       linewidth=0.5, 
                       aes(fill=group), 
                       alpha=0.3, 
                       show.legend = T) +
  scale_fill_manual(values = c("#f2ccac","#a1d5b9","#e1abbc")) # 这里添加或修改颜色
# p2

#install.packages("ggExtra") # 如果尚未安装
# “density”, “histogram”, “boxplot”, “violin”, “densigram”
# 添加边缘箱线图
# 指定图形的宽度、高度和分辨率
pdf(file = "RHYXM_16S.anosim.NMDS.pdf", width=8, height=8) 
## pdf格式
# pdf('marginal.pdf',width = 5,height = 5)
ggMarginal(
  p2,
  type =c('boxplot'),
  margins = 'both',
  size = 3.5,
  groupColour = F,
  groupFill = T
) 
dev.off()

#
# 对于NMDS二维分析，通常认为stress<0.2时有一定的解释意义；当stress<0.1时，可认为是一个好的排序；当 stress<0.05时，则具有很好的代表性。
# stress=0.135

# Adonis,P>0.05,表示不同组的样品之间不存在显著差异
