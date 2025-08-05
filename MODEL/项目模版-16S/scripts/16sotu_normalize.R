rm(list=ls())

library(reshape2)
library(ggplot2)
library(RColorBrewer)
# install.packages('egg')
library(ggalluvial)
library(ggh4x)
library(dplyr)
library(egg)
library(tidyr)

data_raw <- read.csv(
  file = "./otu_taxa_table.xls",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t",
  na.strings = "",
  fill=TRUE
)
df <- data_raw

####获取tax####
df_split <- separate(df, col = taxonomy, into = paste0("col", 1:7), sep = ";")  
rownames(df_split) <- df_split$`#OTU ID`
df_split <- df_split[,-1]
last_7_cols <- df_split[, -(1:(ncol(df_split) - 7))] 
colnames(last_7_cols) <- c('Kingdom','Phylum',
                           'Class','Order',
                           'Family','Genus',
                           'Species')

tax <- last_7_cols



####获取otu####
otu <- data_raw
original_colnames <- colnames(otu) 
cleaned_colnames <- gsub(".bracken.S", "", original_colnames)  
# 将修改后的列名赋值回数据框  
colnames(otu) <- cleaned_colnames 
otu <- otu %>% select(-taxonomy)
rownames(otu) <- otu$`#OTU ID`
otu <- otu[,-1]

library(phyloseq)

otu_table_phy <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
taxonomy_table_phy <- tax_table(as.matrix(tax))    
#合并
physeq <- phyloseq(otu_table_phy, taxonomy_table_phy) 

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
Genus_abundance_with_taxonomy <- Genus_abundance_with_taxonomy[!grepl("g__uncultured", Genus_abundance_with_taxonomy$taxonomy), ] 

# df_with_new_col <- data.frame(new_column = c(rownames(Genus_abundance_with_taxonomy)), df)  
# 保存结果到文件
write.csv(Genus_abundance_with_taxonomy, file = "Genus_abundance_with_taxonomy.csv", row.names = FALSE)
write.table(Genus_abundance_with_taxonomy,file = 'Genus.tmp',sep = "\t", row.names = FALSE, quote = FALSE, na = "")
