#devtools::install_github("jaspershen/metPath")
##########################HDML
setwd('E:\\boyu\\2024-6-30肝脏损伤-李丛舒\\202486丛舒代谢组\\')
library(metPath)
library(ggplot2)
data("hmdb_pathway", package = "metPath")
data("pathbank_pathway", package = "metPath")
hmdb_pathway

P_sig<- read.csv('P_group2 v.s. group1_sig.csv',row.names = 'X')
head(P_sig)

P_ID<- openxlsx::read.xlsx('肝脏-POS-差异离子鉴定(1).xlsx') 
rownames(P_ID)<- P_ID[,1]
head(P_ID)
write.table(P_ID[,4][!is.na(P_ID[,4])],file='name.txt',
            quote = F,row.names = F,col.names = F)

#write.csv(P_sig,file = 'P_sig.csv')
table(rownames(P_sig) %in% rownames(P_ID))
int<- intersect(rownames(P_sig),P_ID[,1])
P_sig<- P_sig[int,]
P_ID<- P_ID[int,]
table(rownames(P_sig)==P_ID[,1])
P_sig$gene<- P_ID[,3]
table(is.na(P_sig$gene))
P_sig<- P_sig[!is.na(P_sig$gene),]
table(P_sig$change)
P_ID<- P_ID[rownames(P_sig),]

d1<- P_sig[grepl('HMDB',P_sig$gene),]
d2<- P_sig[grepl('CSID',P_sig$gene),]
for (i in 1:dim(d1)[1]) {
  d1$gene[i]<- paste0('hmdb:',d1$gene[i])
}
for (i in 1:dim(d2)[1]) {
  temp<- d2$gene[i]
  temp<- strsplit(temp,split = 'CSID')[[1]][2]
  d2$gene[i]<- paste0('chemspider:',temp)
}

write.table(c(d1[d1$change=='Up',]$gene,d2[d2$change=='Up',]$gene),
            file = 'C:\\Users\\wmy\\Desktop\\RaMP_P_up.txt',
            quote = F,col.names = F,row.names = F)
get_pathway_class(hmdb_pathway)
cc<- unique(N_sig[N_sig$change=='Down','gene'])
cc<- c(cc_N_down,cc_P_down)
cc<- cc[grepl('HMDB',cc)]
result = 
  enrich_hmdb(query_id = unique(cc), 
              query_type = "compound", 
              id_type = "HMDB",
              pathway_database = hmdb_pathway,
              only_primary_pathway = TRUE,
              p_cutoff = 1, 
              p_adjust_method = "BH", 
              threads = 3)
result
result<- result@result
result<- result[order(result$p_value),]
top10<- result[1:1,]
top10$text_x <- rep(0.03,1)
top10<- top10[!top10$mapped_number==0,]
top10<- top10[top10$p_value<0.05,]
top10<- top10[order(top10$mapped_percentage),]
nam<- top10$pathway_name
top10$pathway_name<- factor(top10$pathway_name,levels = nam)
write.csv(result[result$p_value<0.05,],file = 'HMDB_down_dayu05.csv')

top10<- read.csv('C:\\Users\\wmy\\Desktop\\mbrole3_disease_up(1).csv')
head(top10)
top10$text_x <- rep(0.03,1)
head(top10)
p2 <- ggplot(data = top10,
             aes(x = mapped_percentage, y = pathway_name)) +
  geom_bar(aes(fill = p_value), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(x = "Ratio",title = "HMDB_down") +
  geom_text(aes(x = text_x, #用新增的重复数组控制文本标签起始位置
                label = pathway_name),
            hjust = 0)+ #hjust=0，左对齐
  theme_classic() 
p2
pdf('HMDB_down_dayu05.pdf',6,3)
dev.off()

###########KEGG
library(metPath)
library(tidyverse)
data("kegg_hsa_pathway", package = "metPath")
kegg_hsa_pathway
remain_idx =
  kegg_hsa_pathway@pathway_class %>%
  unlist() %>%
  stringr::str_detect("Disease") %>%
  `!`() %>%
  which()

pathway_database =
  filter_pathway(object = kegg_hsa_pathway, remain_idx = remain_idx)


P_sig<- read.csv('P_group2 v.s. group1_sig.csv',row.names = 'X')
head(P_sig)

P_ID<- openxlsx::read.xlsx('肝脏-POS-差异离子鉴定(1).xlsx') 
rownames(P_ID)<- P_ID[,1]
head(P_ID)

table(rownames(P_sig) %in% rownames(P_ID))
int<- intersect(rownames(P_sig),P_ID[,1])
P_sig<- P_sig[int,]
P_ID<- P_ID[int,]
table(rownames(P_sig)==P_ID[,1])
P_sig$gene<- P_ID[,3]
table(is.na(P_sig$gene))
P_sig<- P_sig[!is.na(P_sig$gene),]
table(P_sig$change)
P_ID<- P_ID[rownames(P_sig),]

up<- P_ID[rownames(P_sig[P_sig$change=='Up',]),4]
down<- P_ID[rownames(P_sig[P_sig$change=='Down',]),4]

KEGG_ID<- read.csv('P_KEGG.csv')
KEGG_ID<- KEGG_ID[!is.na(KEGG_ID$KEGG),]
up_P<- KEGG_ID[which(KEGG_ID$Query %in% up),]$KEGG
down_P<- KEGG_ID[which(KEGG_ID$Query %in% down),]$KEGG

down<- unique(c(up_N,up_P))

cc<- down
result = 
  enrich_kegg(query_id =unique(cc) , 
              query_type = "compound", 
              id_type = "KEGG",
              pathway_database = pathway_database, 
              p_cutoff = 1, 
              p_adjust_method = "BH", 
              threads = 3)
result
result<- result@result
enrich<- result[order(result$p_value),]
top10<- enrich[1:3,]
write.csv(result[result$p_value<0.05,],file = 'KEGG_up_dayu05.csv')
###ggplot
top10$text_x <- rep(0.03,3)
top10<- top10[!top10$mapped_number==0,]
top10<- top10[top10$p_value<0.05,]
top10<- top10[order(top10$mapped_percentage),]
nam<- top10$pathway_name
top10$pathway_name<- factor(top10$pathway_name,levels = nam)
p2 <- ggplot(data = top10,
             aes(x = mapped_percentage, y = pathway_name)) +
  geom_bar(aes(fill = p_value), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(x = "Ratio", y = "pathway", title = "KEGG_up") +
  geom_text(aes(x = text_x, #用新增的重复数组控制文本标签起始位置
                label = pathway_name),
            hjust = 0)+ #hjust=0，左对齐
  theme_classic() 
p2

pdf('KEGG_updayu05.pdf',8,6)
dev.off()

####################
pa<- read.table('C:\\Users\\wmy\\Desktop\\mbrole3-direct_24032107 (1).tsv',fill = T)

library(ggplot2)
colnames(pa)[8]<- 'In.set'
p<- ggplot(pa,aes(x=Database,y=Annotation))
p+geom_point() #绘制散点图
p+geom_point(aes(size=In.set,color=FDR))+scale_color_gradient(low="red",high="green")+ theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
