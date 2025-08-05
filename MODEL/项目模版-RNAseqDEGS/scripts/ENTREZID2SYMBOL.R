#!/usr/bin/Rscript
####load packages#####
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)

# cmd:Rscript ENTREZID2SYMBOL.R GSE188774 Hs keggentrezid.txt
# args[1] 项目名称：GSE188774
# args[2] 数据库缩写：Hs
# args[3] 包含ENTREZID的文件：keggentrezid.txt

# 查看操作路径
getwd()
# 获取参数
args <- commandArgs(trailingOnly = TRUE)
print(args)
projectname <- args[1]#GSE188774
# 数据库为Hs/Mm
db_short <- args[2]#Hs
if (db_short=='Hs'){
  db='org.Hs.eg.db'
}else if (db_short=='Mm'){
  db='org.Mm.eg.db'
}
#读取处理后的Genes文本文件
Genesfile <- args[3] # 'test/GeneUP-unique.txt'
Genesdf <- read.csv(Genesfile,col.names = 'ENTREZID',header = FALSE)
length(Genesdf$ENTREZID)
Geneslist <- Genesdf$ENTREZID

gs<- bitr(Geneslist, fromType = "ENTREZID", toType="SYMBOL",OrgDb = db)
gs <- gs[!is.na(gs$SYMBOL),]
SYMBOLdf <- as.data.frame(gs$SYMBOL)

write.table(SYMBOLdf, file = sprintf("%s-SYMBOL.txt",projectname ),
            sep = "\t",         # 分隔符
            row.names = FALSE,  # 不寫入行名
            col.names = FALSE,   # 不保留列名
            quote = FALSE,      # 不添加引號
            na = "NA"           # 缺失值標記（自定義為 NA）
)






