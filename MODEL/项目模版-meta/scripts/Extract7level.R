rm(list=ls())

# 加载所需的包
library(reshape2)
library(ggplot2)
library(RColorBrewer)
# install.packages('egg')
library(ggalluvial)
library(ggh4x)
library(dplyr)
library(egg)
library(tidyr)

library(readxl) 
# 加载种水平物种丰度表，并设置列名和分隔符
data_raw <- read.csv(
  file = "./0828.S.count.tsv.all.tmp",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t",
  na.strings = "",
  fill=TRUE
)

# taxonomy <- sapply(strsplit(data$taxonomy, "; g__"), function(x) x[length(x)])
# rownames(data) <- taxonomy
# # 删除第一列OTU
# data <- data[,-1]
# # 删除最后一列taxonomy
# data <- data %>% select(-taxonomy)
# 
# original_colnames <- colnames(data) 
# cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# # 将修改后的列名赋值回数据框  
# colnames(data) <- cleaned_colnames  
# 
# samples <- colnames(data)

otu <- data_raw
original_colnames <- colnames(otu) 
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
colnames(otu) <- cleaned_colnames 
otu <- otu %>% select(-taxonomy)
rownames(otu) <- otu$`#OTU ID`
otu <- otu[,-1]

# 使用separate()函数分割combined_data列，并生成6列新数据  
df_split <- separate(data_raw, col = taxonomy, into = paste0("col", 1:7), sep = ";")  
rownames(df_split) <- df_split$`#OTU ID`
df_split <- df_split[,-1]
last_7_cols <- df_split[, -(1:(ncol(df_split) - 7))] 
colnames(last_7_cols) <- c('Kingdom','Phylum',
                           'Class','Order',
                           'Family','Genus',
                           'Species')

tax <- last_7_cols
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")  



library(phyloseq)
# otu_table <-read.csv("otu_table.csv",header = T,row.names = 1)
# 
# taxonomy_table <- read.csv("taxonomy_table.csv",header = T,row.names = 1)
otu_table_phy <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
taxonomy_table_phy <- tax_table(as.matrix(tax))    
#合并
physeq <- phyloseq(otu_table_phy, taxonomy_table_phy) 

######提取Kingdom#####
Kingdom_abundance <- tax_glom(physeq, taxrank = "Kingdom")
# 提取门的丰度表
Kingdom_abundance_table <- otu_table(Kingdom_abundance)
# 提取门的分类信息
Kingdom_taxonomy <- tax_table(Kingdom_abundance)
# 将丰度表和分类信息结合
Kingdom_abundance_with_taxonomy <- cbind(as.data.frame(Kingdom_abundance_table), as.data.frame(Kingdom_taxonomy))


taxonomy <- Kingdom_abundance_with_taxonomy$Kingdom
# 首先获取df的列数  
ncol_df <- ncol(Kingdom_abundance_with_taxonomy)  
# 然后选择除了最后7列之外的所有列  
Kingdom_abundance_with_taxonomy <- Kingdom_abundance_with_taxonomy[, 1:(ncol_df - 7)]
# df_abun <- t(Kingdom_abundance_with_taxonomy)
Kingdom_abundance_with_taxonomy$taxonomy <- taxonomy

# df_with_new_col <- data.frame(new_column = c(rownames(Kingdom_abundance_with_taxonomy)), df)  
# 保存结果到文件
write.csv(Kingdom_abundance_with_taxonomy, file = "Kingdom_abundance_with_taxonomy.csv")
# data <- Kingdom_abundance_with_taxonomy


data_Kingdom<- read.csv(
  file = "./Kingdom_abundance_with_taxonomy.csv",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = ",",
  na.strings = "",
  fill=TRUE
)
colnames(data_Kingdom)[1] <- "#OTU ID"
write.csv(data_Kingdom, file = "Kingdom_abundance.csv",row.names = FALSE)


######提取Phylum#####
Phylum_abundance <- tax_glom(physeq, taxrank = "Phylum")
# 提取门的丰度表
Phylum_abundance_table <- otu_table(Phylum_abundance)
# 提取门的分类信息
Phylum_taxonomy <- tax_table(Phylum_abundance)
# 将丰度表和分类信息结合
Phylum_abundance_with_taxonomy <- cbind(as.data.frame(Phylum_abundance_table), as.data.frame(Phylum_taxonomy))


taxonomy <- Phylum_abundance_with_taxonomy$Phylum
# 首先获取df的列数  
ncol_df <- ncol(Phylum_abundance_with_taxonomy)  
# 然后选择除了最后7列之外的所有列  
Phylum_abundance_with_taxonomy <- Phylum_abundance_with_taxonomy[, 1:(ncol_df - 7)]
# df_abun <- t(Phylum_abundance_with_taxonomy)
Phylum_abundance_with_taxonomy$taxonomy <- taxonomy

# df_with_new_col <- data.frame(new_column = c(rownames(Phylum_abundance_with_taxonomy)), df)  
# 保存结果到文件
write.csv(Phylum_abundance_with_taxonomy, file = "Phylum_abundance_with_taxonomy.csv")
# data <- Phylum_abundance_with_taxonomy


data_Phylum<- read.csv(
  file = "./Phylum_abundance_with_taxonomy.csv",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = ",",
  na.strings = "",
  fill=TRUE
)
colnames(data_Phylum)[1] <- "#OTU ID"
write.csv(data_Phylum, file = "Phylum_abundance.csv",row.names = FALSE)



######提取Class#####
Class_abundance <- tax_glom(physeq, taxrank = "Class")
# 提取门的丰度表
Class_abundance_table <- otu_table(Class_abundance)
# 提取门的分类信息
Class_taxonomy <- tax_table(Class_abundance)
# 将丰度表和分类信息结合
Class_abundance_with_taxonomy <- cbind(as.data.frame(Class_abundance_table), as.data.frame(Class_taxonomy))


taxonomy <- Class_abundance_with_taxonomy$Class
# 首先获取df的列数  
ncol_df <- ncol(Class_abundance_with_taxonomy)  
# 然后选择除了最后7列之外的所有列  
Class_abundance_with_taxonomy <- Class_abundance_with_taxonomy[, 1:(ncol_df - 7)]
# df_abun <- t(Class_abundance_with_taxonomy)
Class_abundance_with_taxonomy$taxonomy <- taxonomy

# df_with_new_col <- data.frame(new_column = c(rownames(Class_abundance_with_taxonomy)), df)  
# 保存结果到文件
write.csv(Class_abundance_with_taxonomy, file = "Class_abundance_with_taxonomy.csv")
# data <- Class_abundance_with_taxonomy


data_Class<- read.csv(
  file = "./Class_abundance_with_taxonomy.csv",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = ",",
  na.strings = "",
  fill=TRUE
)
colnames(data_Class)[1] <- "#OTU ID"
write.csv(data_Class, file = "Class_abundance.csv",row.names = FALSE)

######提取Order#####
Order_abundance <- tax_glom(physeq, taxrank = "Order")
# 提取门的丰度表
Order_abundance_table <- otu_table(Order_abundance)
# 提取门的分类信息
Order_taxonomy <- tax_table(Order_abundance)
# 将丰度表和分类信息结合
Order_abundance_with_taxonomy <- cbind(as.data.frame(Order_abundance_table), as.data.frame(Order_taxonomy))


taxonomy <- Order_abundance_with_taxonomy$Order
# 首先获取df的列数  
ncol_df <- ncol(Order_abundance_with_taxonomy)  
# 然后选择除了最后7列之外的所有列  
Order_abundance_with_taxonomy <- Order_abundance_with_taxonomy[, 1:(ncol_df - 7)]
# df_abun <- t(Order_abundance_with_taxonomy)
Order_abundance_with_taxonomy$taxonomy <- taxonomy

# df_with_new_col <- data.frame(new_column = c(rownames(Order_abundance_with_taxonomy)), df)  
# 保存结果到文件
write.csv(Order_abundance_with_taxonomy, file = "Order_abundance_with_taxonomy.csv")
# data <- Order_abundance_with_taxonomy


data_Order<- read.csv(
  file = "./Order_abundance_with_taxonomy.csv",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = ",",
  na.strings = "",
  fill=TRUE
)
colnames(data_Order)[1] <- "#OTU ID"
write.csv(data_Order, file = "Order_abundance.csv",row.names = FALSE)
######提取Family#####
Family_abundance <- tax_glom(physeq, taxrank = "Family")
# 提取门的丰度表
Family_abundance_table <- otu_table(Family_abundance)
# 提取门的分类信息
Family_taxonomy <- tax_table(Family_abundance)
# 将丰度表和分类信息结合
Family_abundance_with_taxonomy <- cbind(as.data.frame(Family_abundance_table), as.data.frame(Family_taxonomy))


taxonomy <- Family_abundance_with_taxonomy$Family
# 首先获取df的列数  
ncol_df <- ncol(Family_abundance_with_taxonomy)  
# 然后选择除了最后7列之外的所有列  
Family_abundance_with_taxonomy <- Family_abundance_with_taxonomy[, 1:(ncol_df - 7)]
# df_abun <- t(Family_abundance_with_taxonomy)
Family_abundance_with_taxonomy$taxonomy <- taxonomy

# df_with_new_col <- data.frame(new_column = c(rownames(Family_abundance_with_taxonomy)), df)  
# 保存结果到文件
write.csv(Family_abundance_with_taxonomy, file = "Family_abundance_with_taxonomy.csv")
# data <- Family_abundance_with_taxonomy


data_Family<- read.csv(
  file = "./Family_abundance_with_taxonomy.csv",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = ",",
  na.strings = "",
  fill=TRUE
)
colnames(data_Family)[1] <- "#OTU ID"
write.csv(data_Family, file = "Family_abundance.csv",row.names = FALSE)

######提取Genus#####
Genus_abundance <- tax_glom(physeq, taxrank = "Genus")
# 提取门的丰度表
Genus_abundance_table <- otu_table(Genus_abundance)
# 提取门的分类信息
Genus_taxonomy <- tax_table(Genus_abundance)
# 将丰度表和分类信息结合
Genus_abundance_with_taxonomy <- cbind(as.data.frame(Genus_abundance_table), as.data.frame(Genus_taxonomy))


taxonomy <- Genus_abundance_with_taxonomy$Genus
# 首先获取df的列数  
ncol_df <- ncol(Genus_abundance_with_taxonomy)  
# 然后选择除了最后7列之外的所有列  
Genus_abundance_with_taxonomy <- Genus_abundance_with_taxonomy[, 1:(ncol_df - 7)]
# df_abun <- t(Genus_abundance_with_taxonomy)
Genus_abundance_with_taxonomy$taxonomy <- taxonomy

# df_with_new_col <- data.frame(new_column = c(rownames(Genus_abundance_with_taxonomy)), df)  
# 保存结果到文件
write.csv(Genus_abundance_with_taxonomy, file = "Genus_abundance_with_taxonomy.csv")
# data <- Genus_abundance_with_taxonomy


data_Genus<- read.csv(
  file = "./Genus_abundance_with_taxonomy.csv",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = ",",
  na.strings = "",
  fill=TRUE
)
colnames(data_Genus)[1] <- "#OTU ID"
write.csv(data_Genus, file = "Genus_abundance.csv",row.names = FALSE)

######提取Species#####
Species_abundance <- tax_glom(physeq, taxrank = "Species")
# 提取门的丰度表
Species_abundance_table <- otu_table(Species_abundance)
# 提取门的分类信息
Species_taxonomy <- tax_table(Species_abundance)
# 将丰度表和分类信息结合
Species_abundance_with_taxonomy <- cbind(as.data.frame(Species_abundance_table), as.data.frame(Species_taxonomy))


taxonomy <- Species_abundance_with_taxonomy$Species
# rownames(Species_abundance_with_taxonomy) <- taxonomy
# 首先获取df的列数  
ncol_df <- ncol(Species_abundance_with_taxonomy)  
# 然后选择除了最后7列之外的所有列  
Species_abundance_with_taxonomy <- Species_abundance_with_taxonomy[, 1:(ncol_df - 7)]
# df_abun <- t(Species_abundance_with_taxonomy)
Species_abundance_with_taxonomy$taxonomy <- taxonomy

# df_with_new_col <- data.frame(new_column = c(rownames(Species_abundance_with_taxonomy)), df)  
# 保存结果到文件
write.csv(Species_abundance_with_taxonomy, file = "Species_abundance_with_taxonomy.csv")
# data <- Species_abundance_with_taxonomy


data_Species<- read.csv(
  file = "./Species_abundance_with_taxonomy.csv",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = ",",
  na.strings = "",
  fill=TRUE
)
colnames(data_Species)[1] <- "#OTU ID"
write.csv(data_Species, file = "Species_abundance.csv",row.names = FALSE)



####新的提取方法####
rare = read.table("rare.txt", header=T, sep="\t", row.names=1, comment.char="")

# 提取phylum，删除无用
phylum=c()
delete=c()
count = 1
for(i in rare[, length(rare)])
{
  tmp = unlist(strsplit(as.character(i), split="; |__"))
  phylum = c(phylum, tmp[4])
  if(is.na(tmp[4]) | tmp[4] == "")
  {
    delete = c(delete, count)
  }
  count = count + 1
}

rare$phylum = phylum
rare = rare[, -(ncol(rare)-1)]  # 删除注释
rare2 = rare[-delete, ]  # 删除NA和""
rare2 = rare2[order(rare2$phylum, decreasing=F),]
plist = unique(rare2$phylum)

# 合并行，相同门
rare3 = data.frame(apply(rare2[rare2$phylum==plist[1], c(1:(ncol(rare2)-1))], 2, sum))
colnames(rare3)[1] = plist[1]
for(i in 2:length(plist))
{
  tmp = apply(rare2[rare2$phylum==plist[i], c(1:(ncol(rare2)-1))], 2, sum)
  rare3 = cbind(rare3, tmp)
  colnames(rare3)[i] = plist[i]
}
rare3 = data.frame(t(rare3))

# 数据归一化
norm = rare3
sample_sum = apply(rare3, 2, sum)
for(i in 1:nrow(rare3))
{
  for(j in 1:ncol(rare3))
  {
    norm[i, j] = rare3[i, j]/sample_sum[j]
  }
}
apply(norm, 2, sum)  # 检测
write.csv(norm, file="phylum_rare_norm.csv")


