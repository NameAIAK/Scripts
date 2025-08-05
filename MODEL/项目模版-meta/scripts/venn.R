rm(list=ls())
library(tidyverse)
# input OTUs 
data  <-  read.delim("./0828.S.count.tsv.all.tmp", header=T, sep="\t", row.names=1, stringsAsFactors = FALSE)

metadata <- read.delim("./0828.group.txt", header=T, sep="\t", stringsAsFactors = FALSE)
tax <- data.frame(data[,'taxonomy'])
rownames(tax) <- rownames(data)
colnames(tax) <- c('taxonomy')

otu_table <- data %>% select(-taxonomy)
original_colnames <- colnames(otu_table) 
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
colnames(otu_table) <- cleaned_colnames 

colnames(metadata) <- c('Sample','Group')
# check data
index <- metadata$Sample %in% colnames(otu_table)
metadata <- metadata[index, ]
otu_table <- otu_table[ ,colnames(otu_table)] 

# matrix transpose
transp_otu <- as.data.frame(t(otu_table))
transp_otu$Sample <- rownames(transp_otu)

transp_otu <- merge(metadata, transp_otu, by="Sample")
# transp_otu[c(1:3), c(1:3)]
transp_otu <- transp_otu[,c(-1)]
# transp_otu[c(1:3), c(1:3)]
# install.packages('tidyverse')

transp_otu %>% group_by(Group) %>% 
  summarise_all(mean) -> group_otu

# threshold value
## VennDiagram 要求输入每个圆圈的元素名列表，可以先进行筛选然后再可视化。
## 当然，也可以在可视化的时候进行筛选。
venn_otu <- as.data.frame(group_otu)
rownames(venn_otu) <- venn_otu$Group
venn_otu <-as.data.frame(t(venn_otu[,-1]))
# backups 
venn_otu2 <- venn_otu

venn_otu[venn_otu > 0] <- 1

venn_otu$id <- rownames(venn_otu) 
tax$id <- rownames(tax)

merge <- merge(venn_otu, tax, by="id")
rownames(merge) <- merge$taxonomy
merge <- merge[,2:3]




# visualization
# install.packages("VennDiagram")
library(VennDiagram)

# p1 <- venn.diagram(list(CONTROL=row.names(venn_otu[venn_otu$CONTROL==1, ]),
#                            OSAHS=row.names(venn_otu[venn_otu$OSAHS==1, ])),
#                       fill = c("dodgerblue", "goldenrod1"),
#                       alpha = 0.5,cex = 2,
#                       cat.cex = 2,cat.fontface = "bold",lwd=3,
#                       resolution = 300,filename = NULL)
# 
# pdf('0828.venn.pdf', width = 10, height = 10)
# grid.draw(p1)
# dev.off()


# p2 <- venn.diagram(
#   x=list(C=row.names(venn_otu2[venn_otu2$C>0, ]),
#          O=row.names(venn_otu2[venn_otu2$O>0, ])),
#   filename = "venn2.png", lwd = 3, alpha = 0.6,
#   label.col = "white", cex = 1.5,
#   fill = c("#eb507e", "#2f90b9"), 
#   cat.col = c("#eb507e", "#2f90b9"),
#   fontfamily = "serif", fontface = "bold",
#   cat.fontfamily = "serif",cat.fontface = "bold",
#   margin = 0.05)



venn_list <- list(OSAHS = row.names(merge[merge$OSAHS==1, ]), CONTROL = row.names(merge[merge$CONTROL==1, ]))

p_venn<- venn.diagram(venn_list, imagetype = 'png', 
                            fill = c("dodgerblue", "goldenrod1"), alpha = 0.50, cat.col = rep('black', 2), 
                            col = 'black', cex = 1.5, fontfamily = 'serif', 
                            cat.cex = 1.5, cat.fontfamily = 'serif',filename = NULL)

pdf('species_venn.pdf', width = 10, height = 10)
grid.draw(p_venn)
dev.off()

inter <- get.venn.partitions(venn_list)

control <- data.frame(inter[[2,'..values..']])
colnames(control) <- c('taxonomy')
control_data <- merge(data,control,by='taxonomy')
rownames(control_data) <- control_data$taxonomy
control_tax <- control_data$taxonomy
control_data <- control_data %>% select(-taxonomy)
control_data$sum <- rowSums(control_data)
# colSums(df_normalized)
# 按每行的和降序排列
data1 <- control_data[order(control_data$sum, decreasing = TRUE), ]
write.table(data1, 'control_venn_inter.csv', row.names = TRUE, col.names = TRUE,sep = ',', quote = FALSE)

rownames(data1)[1:10]


osahs <- data.frame(inter[[3,'..values..']])
colnames(osahs) <- c('taxonomy')
osahs_data <- merge(data,osahs,by='taxonomy')
rownames(osahs_data) <- osahs_data$taxonomy
osahs_tax <- osahs_data$taxonomy
osahs_data <- osahs_data %>% select(-taxonomy)
osahs_data$sum <- rowSums(osahs_data)
# colSums(df_normalized)
# 按每行的和降序排列
data2 <- osahs_data[order(osahs_data$sum, decreasing = TRUE), ]
write.table(data2, 'osahs_venn_inter.csv', row.names = TRUE, col.names = TRUE,sep = ',', quote = FALSE)

rownames(data2)[1:10]



# for (i in 1:nrow(inter)) inter2[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
# write.table(inter2[-c(3, 4)], 'species_venn_inter.txt', row.names = FALSE, sep = '\t', quote = FALSE)

