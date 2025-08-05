# First
datadir <- './data1'
files <- list.files(datadir)

df <- data.frame(Compound = character(),Accepted.Compound.ID = character(),Accepted.Description=character())
for (f in files){
  filepath <- file.path(datadir,f)
  data <- read.csv(filepath)
  df <- rbind(df,data)
  df <- unique(df)
}
write.csv(df,'MetaboliteName20250425.csv',row.names = FALSE)

library(dplyr)
#add
df <- read.csv('./MetaboliteName20250425.csv')
newdatadir <- './data2'
newfiles <- list.files(newdatadir)
for (newf in newfiles){
  filepath <- file.path(newdatadir,newf)
  data <- read.csv(filepath)
  df <- rbind(df,data)
  df <- unique(df)
}

# 统计每个 Compound 的频率
compound_counts <- table(df$Compound)
# 提取出现次数 >1 的 Compound（即重复项）
repeated_compounds <- names(compound_counts[compound_counts > 1])

df <- df %>%
  filter(!is.na(Accepted.Description) & df$Accepted.Description != "")
df <- df[!duplicated(df$Compound), ]
write.csv(df,'MetaboliteName20250425V2.csv',row.names = FALSE)

