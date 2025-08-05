library(ggplot2)

d1<- read.table('RHvsM_2.txt')
colnames(d1)
d1$gene<- rownames(d1)

# 设置p_value和logFC的阈值
cut_off_FDR = 0.05  #统计显著性
cut_off_logFC = 0.5           #差异倍数值
d1$change<- 'Stable'
# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
d1$change = ifelse(d1$FDR< cut_off_FDR & abs(d1$logFC) > cut_off_logFC, 
                        ifelse(d1$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(d1)
table(d1$change)

p <- ggplot(
  # 数据、映射、颜色
  d1, aes(x = logFC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_FDR),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  scale_x_continuous(limits = c(-5, 5))+
  scale_y_continuous(limits = c(0,30))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('RH v.s. M')
p

# 将需要标记的基因放置在label列(logFC >= 5)
library(ggrepel)
d1$label <- ifelse(d1$PValue < cut_off_FDR & abs(d1$logFC) >= 0.3,
                        as.character(d1$gene), "")


p + geom_label_repel(data = d1, aes(x = d1$logFC, 
                                         y = -log10(d1$PValue), 
                                         label = label),
                     size = 3, box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"), 
                     segment.color = "black", 
                     show.legend = FALSE)

library(ggvenn)
table(MC$change)
up<- list(
  'M vs C'=MC[MC$change=='Up',]$gene,
  'RH vs C'=RHC[RHC$change=='Up',]$gene,
  'RH vs M'=RHM[RHM$change=='Up',]$gene
)
down<- list(
  'M vs C'=MC[MC$change=='Down',]$gene,
  'RH vs C'=RHC[RHC$change=='Down',]$gene,
  'RH vs M'=RHM[RHM$change=='Down',]$gene
)
#up<- list('M vs C'=MC[MC$change=='Up',]$gene,'M vs B'=MB[MB$change=='Up',]$gene)
#down<- list('M vs C'=MC[MC$change=='Down',]$gene,'M vs B'=MB[MB$change=='Down',]$gene)

ggvenn(up,c('M vs C', 'M vs B'),show_percentage = T,
       stroke_color = "white",
       fill_color = c("#BC3C29FF","#0072B5FF","#E18727FF","#FFDC91FF",
                      "#EE4C97FF","#20854EFF","#7876B1FF","#6F99ADFF" ),
       set_name_color =c("#E41A1C","#1E90FF",'#20854EFF'))


