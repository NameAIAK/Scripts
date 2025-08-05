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


####16S属数据合并####

data_raw <- read.csv(
  file = "./ASV_table_even_tax.xls",
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
parts <- strsplit(df$taxonomy, "\\|")
num_parts <- sapply(parts, length)

# 过滤掉分割后少于6个部分的行
df_filtered <- df[num_parts >= 6, ]

df_filtered_new <- df_filtered
# 对过滤后的数据框，保留前6个部分并用;连接
df_filtered_new$new_taxonomy <- sapply(strsplit(df_filtered$taxonomy, "\\|"), function(x) paste0(head(x, 5), collapse = ";"))
df_filtered_new[1,'taxonomy']
df_filtered_new <- df_filtered_new%>% select(-taxonomy)
colnames(df_filtered_new)[32] <- 'taxonomy'
ASV <- df_filtered_new$ASV_ID
df_filtered_new <- df_filtered_new[,-1]

grouped_df_all <- df_filtered_new %>%
  group_by(taxonomy) %>%
  summarise_all(sum)


write.table(grouped_df_all, file = "grouped_df_all.F.count.tsv.all.tmp", sep = "\t", row.names = FALSE, quote = FALSE, na = "")





