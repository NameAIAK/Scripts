rm(list=ls())
####KO差异分析####
# 设置p值和log2FoldChange阈值
pvalue <- 0.05
log2fc <- 0


# 加载DESeq2和EnhancedVolcano包
library("DESeq2")
library("EnhancedVolcano")
library('dplyr')  
# install.packages('dplyr')
# setwd('./module/')
count_tmp <- read.csv(
  file = "./0828.KO_samples.xls",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t",
  na.strings = "",
  fill=TRUE
)

# KO_name <- sapply(strsplit(count_tmp$KO_name, "; g__"), function(x) x[length(x)])
class(count_tmp)
# install.packages("data.table")
library(data.table)  
setDT(count_tmp) # 将df转换为data.table格式  
# 使用unique  
df_unique <- unique(count_tmp, by = "KO", fromLast = FALSE)  
count_tmp <- as.data.frame(df_unique)
count_tmp <- na.omit(count_tmp)

rownames(count_tmp) <- count_tmp$KO
count_tmp <- count_tmp[,-c(1,2,3)]

library(readxl) 
# id2 <- read_excel("id2.xlsx")
id <- read_excel("0828.id.xlsx")
# colnames(count_tmp) <- id$id
count_tmp <- count_tmp[,c(id$sample)]

map <- read.table(
  file = "0828.group.txt",
  comment.char = "",
  row.names = 1,
  na.strings = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t"
)

# # 从样本信息中提取Category1信息
# category <- map["Category1"]
# # 移除NA值
# category <- na.omit(category)
# rownames(category) <- id$id
# category <- category[c(id2$id),]
# category <- as.data.frame(category)
# rownames(category) <- id2$id
# colnames(category) <- c('Category1')
# # 在计数数据中加1，以避免在后续步骤中对零取对数
# df <- count_tmp
# df <- df + 1
# category$Category1 <- factor(category$Category1)
map2 <- map[id$sample,]
map2 <- as.data.frame(map2)
rownames(map2) <- id$sample
colnames(map2) <- c('Category1')
map2 <- na.omit(map2)
# 在计数数据中加1，以避免在后续步骤中对零取对数
df_abun <- count_tmp
df_abun <- df_abun + 1
map2$Category1 <- factor(map2$Category1)

df_abun <- na.omit(df_abun)
# 从计数数据和样本信息创建DESeqDataSet对象
dds <- DESeqDataSetFromMatrix(
  countData = df_abun,
  colData = map2,
  design = ~Category1 # 比较的分组方案
)

# 将负二项分布模型拟合到数据并进行差异表达测试
fit <- DESeq(dds, quiet = FALSE)
# 提取差异表达分析的结果，将对比设置为将O组与C组进行比较
contrast = c("Category1", "OSAHS", "CONTROL")
df <- results(
  fit,
  contrast = contrast
)

# 将结果转换为数据框架，并添加一个Features列，其中包含行名
df <- data.frame(
  Features = rownames(df),
  df,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 设置要在火山图的y轴上使用的列，padj列包含校正后的p值
pcol <- "padj"
# 按p值和log2FoldChange绝对值对数据框进行排序
df <- df[
  order(
    df[, pcol],
    abs(df$log2FoldChange),
    decreasing = c(FALSE, TRUE)
  ),]

#####提取差异基因#####
pvalue <- 0.05
log2fc <- 0
pcol <- "padj"

df_up <- df[df$log2FoldChange > log2fc & df[, pcol] < pvalue,]
df_up <- df_up[complete.cases(df_up), ] 
df_upf <- df_up$Features
df_upf
df_down <- df[df$log2FoldChange < -log2fc & df[, pcol] < pvalue,]
df_down <- df_down[complete.cases(df_down), ] 
df_downf <- df_down$Features
df_downf
rbind_df <- rbind(df_up, df_down) 
variables <- rbind_df$Features

df1 <- t(count_tmp[c(variables),])


library(clusterProfiler)
#####导入KOID，富集分析####
# KO_list=variables
#进行富集分析
KO_result_up=enrichKEGG(df_upf,
                     organism = "ko", #物种选择ko
                     pvalueCutoff = 0.05,  #p值cutoff
                     pAdjustMethod = "BH", #FDR矫正p值
                     qvalueCutoff = 0.2,   #q值cutoff
)
#保存富集分析结果

KEGG_result_up <- data.frame(KO_result_up@result)
KEGG_result_up <- KEGG_result_up[KEGG_result_up$pvalue<0.05,]
write.csv(file="combine_KEGG_result_up.csv",KEGG_result_up,row.names=F)

KO_result_down=enrichKEGG(df_downf,
                        organism = "ko", #物种选择ko
                        pvalueCutoff = 0.05,  #p值cutoff
                        pAdjustMethod = "BH", #FDR矫正p值
                        qvalueCutoff = 0.2,   #q值cutoff
)
#保存富集分析结果
# write.csv(file="KO_enrichment_results.csv",KO_result@result,row.names=F)
KEGG_result_down <- data.frame(KO_result_down@result)
KEGG_result_down <- KEGG_result_down[KEGG_result_down$pvalue<0.05,]
write.csv(file="combine_KEGG_result_down.csv",KEGG_result_down,row.names=F)
# 加载所需要的R包
# rm(list = ls())
# library(XML)
# library(RCurl)
# library(tidyverse)
# library(ggplot2)
# library(magrittr)
# library(clusterProfiler)
#第一步下载KEGG数据库信息
# extractCompounds <- function(pathwayId) {
#   compoundUrl <- paste0("https://www.genome.jp/dbget-bin/get_linkdb?-t+compound+path:", pathwayId)
#   compoundDoc <- htmlParse(getURL(compoundUrl), encoding = "utf-8")
#   compoundLinks <- getNodeSet(compoundDoc, "/html/body/pre/a")
#   compoundIds <- sapply(compoundLinks, function(node) xmlGetAttr(node, "href"))
#   compoundNames <- sapply(getNodeSet(compoundDoc, "/html/body/pre/text()"), xmlValue)[-1]
#   data.frame(compoundId = paste(compoundIds, collapse = ";"), compoundName = paste(compoundNames, collapse = ";"))
# }
# 
# # Main process
# keggUrl <- "https://www.genome.jp/kegg/pathway.html#global"
# keggDoc <- htmlParse(getURL(keggUrl), encoding = "UTF-8")
# 
# pathwayLinks <- getNodeSet(keggDoc, "//a[@href]")
# pathwayIds <- sapply(pathwayLinks[65:276], function(node) gsub("/pathway/", "", xmlGetAttr(node, "href")))
# pathwayNames <- sapply(pathwayLinks[65:276], xmlValue)
# 
# # Applying extractCompounds function to each pathwayId
# compoundDataList <- Map(extractCompounds, pathwayIds)
# pathwayData <- do.call(rbind, compoundDataList)
# 
# # Combine all data into a single data frame
# finalData <- data.frame(pathwayId = pathwayIds, pathwayName = pathwayNames, pathwayData)
# finalData %>% write_csv(paste0("KeggAllcompounds-",Sys.Date(),".csv"))
# finalData2 <- finalData[finalData$compoundId != "", ]

#' 分别采用for和map的方式将结果进行整理
#' 整理成为result_data和result_data2
#' 采用for循环的方式 
#' result_data <- data.frame()
#' nrow(finalData2)
#' for (i in 1:nrow(finalData2)) {
#'   cid <- finalData2$compoundId[i]
#'   extracted_cid <- str_extract_all(cid, "C\\d+")
#'   CID <- unlist(extracted_cid)
#'   CName <- finalData2$compoundName[i]
#'   split_CName <- strsplit(CName, "\n;")
#'   CompoundName <- lapply(split_CName[[1]], trimws) %>% unlist()
#'   pathway <- cbind(CID, CompoundName)
#'   pathwayId <- rep(finalData2$pathwayId[i], nrow(pathway))
#'   pathwayName <- rep(finalData2$pathwayName[i], nrow(pathway))
#'   dat <- cbind(pathway, pathwayId, pathwayName)
#'   result_data <- rbind(result_data, dat)
#' }
#' 
#' #' 采用map的方式
#' # Define a function to process each row
#' process_row <- function(row) {
#'   cid <- row$compoundId
#'   extracted_cid <- str_extract_all(cid, "C\\d+")
#'   CID <- unlist(extracted_cid)
#'   
#'   CName <- row$compoundName
#'   split_CName <- strsplit(CName, "\n;")
#'   CompoundName <- lapply(split_CName[[1]], trimws) %>% unlist()
#'   
#'   # Ensure the result is a data frame
#'   if (length(CID) > 0 && length(CompoundName) > 0) {
#'     pathway <- data.frame(CID, CompoundName, stringsAsFactors = FALSE)
#'     pathwayId <- rep(row$pathwayId, nrow(pathway))
#'     pathwayName <- rep(row$pathwayName, nrow(pathway))
#'     return(data.frame(pathway, pathwayId, pathwayName, stringsAsFactors = FALSE))
#'   } else {
#'     return(data.frame(CID = character(), CompoundName = character(), pathwayId = row$pathwayId, pathwayName = row$pathwayName, stringsAsFactors = FALSE))
#'   }
#' }
#' 
#' # Apply the process_row function to each row of finalData using map_df
#' result_data2 <- map_df(seq_len(nrow(finalData2)), ~process_row(finalData2[.x, ]))
#' 
#' result_data %>% write_csv(paste0("keggAllCompoundReshapedData2-",Sys.Date(),".csv")) #将整理好的结果进行储存，后面可以直接读取进来用，不用重复跑前面的代码，毕竟也需要时间的
#' result_data2 %>% write_csv(paste0("keggAllCompoundReshapedData22-",Sys.Date(),".csv")

#####metabolite#####

metabolite_row <- read_excel("metabolite_data.xlsx")
ncol_metabolite <- ncol(metabolite_row)
metabolite <- metabolite_row[,c(2,11,17:ncol_metabolite)]
keggannotation <- read.csv('keggannotation.csv')
keggannotation <- keggannotation[,c(2,3)]
# write.csv(keggannotation,file='keggannotation.csv')

allkeggid <- metabolite$KEGG
allkeggid <- data.frame(allkeggid)
colnames(allkeggid) <- c('ID')
total <- right_join(keggannotation,allkeggid,by="ID") %>% select(2,1)
#####拆分pos/neg#####
mb_pos <- metabolite[c(metabolite$posneg=='pos'),]
mb_pos <- mb_pos[,c('KEGG',id$id)]
mb_pos <- na.omit(mb_pos)
rownames_mb_pos <- mb_pos$KEGG
mb_pos <- mb_pos[,-1]
rownames(mb_pos) <- rownames_mb_pos
df_mb_pos <- t(mb_pos)


mb_neg <- metabolite[c(metabolite$posneg=='neg'),]
mb_neg <- mb_neg[,c('KEGG',id$id)]
mb_neg <- na.omit(mb_neg)
rownames_mb_neg <-mb_neg$KEGG
mb_neg <- mb_neg[,-1]
rownames(mb_neg) <- rownames_mb_neg
df_mb_neg <- t(mb_neg)


####pos差异分析####
#if (!requireNamespace("BiocManager", quietly = TRUE))
# BiocManager::install("edgeR")

library(limma)
library(edgeR)

dge <- DGEList(counts=mb_pos)
group.list=map2$Category1
group.list=factor(group.list)
group.list=relevel(group.list,ref = "CONTROL")

design <- model.matrix(~0+group.list)
rownames(design)<-colnames(dge)
colnames(design)<-levels(group.list)
dge <- calcNormFactors(dge)   
v <- voom(dge,design, normalize="quantile")   
fit <- lmFit(v, design)    
constrasts = paste(rev(levels(group.list)),collapse = "-")
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 

fit2=contrasts.fit(fit,cont.matrix)   
fit2=eBayes(fit2)   
DEG = topTable(fit2, coef=constrasts,sort.by = "P", n=Inf)

DEG2=topTable(fit2, coef=1, adjust="BH")  
DEG = na.omit(DEG)   
write.csv(DEG,file = 'combine.P_C v.s. O.csv')

pos_diff_mb <- DEG[DEG[, "P.Value"] < 0.05,]
variables_pos <- rownames(pos_diff_mb)

p_mb_up <- pos_diff_mb[pos_diff_mb$logFC>0,]
variables_pos_up <- rownames(p_mb_up)

p_mb_down <- pos_diff_mb[pos_diff_mb$logFC<0,]
variables_pos_down <- rownames(p_mb_down)


df_pos <- df_mb_pos[,c(variables_pos)]
####neg差异分析#####

library(limma)
library(edgeR)

dge <- DGEList(counts=mb_neg)
group.list=map2$Category1
group.list=factor(group.list)
group.list=relevel(group.list,ref = "CONTROL")

design <- model.matrix(~0+group.list)
rownames(design)<-colnames(dge)
colnames(design)<-levels(group.list)
dge <- calcNormFactors(dge)   
v <- voom(dge,design, normalize="quantile")   
fit <- lmFit(v, design)    
constrasts = paste(rev(levels(group.list)),collapse = "-")
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 

fit2=contrasts.fit(fit,cont.matrix)   
fit2=eBayes(fit2)   
DEG = topTable(fit2, coef=constrasts,sort.by = "P", n=Inf)

DEG2=topTable(fit2, coef=1, adjust="BH")  
DEG = na.omit(DEG)   
write.csv(DEG,file = 'P_C v.s. O.csv')

neg_diff_mb <- DEG[DEG[, "P.Value"] < 0.05,]
variables_neg <- rownames(neg_diff_mb)

n_mb_up <- neg_diff_mb[neg_diff_mb$logFC>0,]
variables_neg_up<- rownames(n_mb_up)

n_mb_down <- neg_diff_mb[neg_diff_mb$logFC<0,]
variables_neg_down <- rownames(n_mb_down)


df_neg <- df_mb_neg[,c(variables_neg)]


# ######gene_pos相关分析#######
# library(psych)  #psych包用于计算相关性、p值等信息
# library(pheatmap) #用于绘制热图
# library(reshape2) #reshape2包用于输出数据的整合处理
# 
# res <- corr.test(df1,df_pos,method = "spearman",alpha = 0.05) #method可选“pearson”、“spearman”、“kendall”
# result_p <- res$p #提取p值
# result_r <- res$r #提取cor值
# p.out<-cbind(rownames(result_p),result_p) #输出p值
# r.out<-cbind(rownames(result_r),result_r) #输出cor值
# write.csv(p.out,"pout.csv",row.names = F)
# write.csv(r.out,"rout.csv",row.names = F)
# 
# df_pr <-melt(result_r,value.name="cor")
# df_pr$pvalue <-as.vector(result_p)  #宽格式转长格式
# df_filt <- subset(df_pr,abs(df_pr$cor)>0.75&df_pr$pvalue<0.05)#筛选
# 
# colnames(df_filt) <- c("from", "to", "value",'p')
# library(circlize)
# ## 简单和弦图绘制
# pdf("KO_metabolite.pdf", width = 8, height = 6)
# chordDiagram(df_filt,annotationTrack = "grid")
# for(i in get.all.sector.index()) {
#   xlim = get.cell.meta.data("xlim", sector.index = i, track.index = 1)
#   circos.text(x = mean(xlim), y = 0.5, sector.index = i, 
#               facing = "clockwise", niceFacing = TRUE, cex = 1, 
#               col = "black", labels = i, track.index = 1)
# }
# circos.clear()
# dev.off()


####共同通路#####

library(XML)
library(RCurl)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(clusterProfiler)
# diffkeggid_up=append(variables_pos_up,variables_neg_up)
# # 富集分析
# mb_result_up <- clusterProfiler::enricher(gene = diffkeggid_up,TERM2GENE = total,minGSSize = 1,pvalueCutoff = 1,qvalueCutoff = 1)
# mb_result_up <- mb_result_up@result
# mb_result_up <- mb_result_up[mb_result_up$pvalue<0.05,]
# # 结果导出
# write.csv(as.data.frame(mb_result_up) %>% select(-1,-2),
#           file = "combine_mb_up_KEGG_enrichment_result.csv")
# 
# 
# diffkeggid_down=append(variables_pos_down,variables_neg_down)
# # 富集分析
# mb_result_down <- clusterProfiler::enricher(gene = diffkeggid_down,TERM2GENE = total,minGSSize = 1,pvalueCutoff = 1,qvalueCutoff = 1)
# mb_result_down <- mb_result_down@result
# mb_result_down <- mb_result_down[mb_result_down$pvalue<0.05,]
# # 结果导出
# write.csv(as.data.frame(mb_result_down) %>% select(-1,-2),
#           file = "combine_mb_down_KEGG_enrichment_result.csv")



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

result_pos_up= 
  enrich_kegg(query_id =unique(variables_pos_up) , 
              query_type = "compound", 
              id_type = "KEGG",
              pathway_database = pathway_database, 
              p_cutoff = 1, 
              p_adjust_method = "BH", 
              threads = 3)
result_pos_up<- result_pos_up@result
result_pos_up <- result_pos_up[result_pos_up$p_value<0.05,]


result_neg_up = 
  enrich_kegg(query_id =unique(variables_neg_up) , 
              query_type = "compound", 
              id_type = "KEGG",
              pathway_database = pathway_database, 
              p_cutoff = 1, 
              p_adjust_method = "BH", 
              threads = 3)
result_neg_up <- result_neg_up@result
result_neg_up <- result_neg_up[result_neg_up$p_value<0.05,]

result_all_up_pathway <- unique(append(result_pos_up$pathway_name,result_neg_up$pathway_name))

result_pos_down= 
  enrich_kegg(query_id =unique(variables_pos_down) , 
              query_type = "compound", 
              id_type = "KEGG",
              pathway_database = pathway_database, 
              p_cutoff = 1, 
              p_adjust_method = "BH", 
              threads = 3)
result_pos_down<- result_pos_down@result
result_pos_down <- result_pos_down[result_pos_down$p_value<0.05,]


result_neg_down = 
  enrich_kegg(query_id =unique(variables_neg_down) , 
              query_type = "compound", 
              id_type = "KEGG",
              pathway_database = pathway_database, 
              p_cutoff = 1, 
              p_adjust_method = "BH", 
              threads = 3)
result_neg_down <- result_neg_down@result
result_neg_down <- result_neg_down[result_neg_down$p_value<0.05,]

result_all_down_pathway <- unique(append(result_pos_down$pathway_name,result_neg_down$pathway_name))
# visualization
# install.packages("VennDiagram")
library(VennDiagram)

venn_list_up <- list(metagenomic = KEGG_result_up$Description, metabolite = result_all_up_pathway)

p_venn_up <- venn.diagram(venn_list_up, imagetype = 'png', 
             fill = c("dodgerblue", "goldenrod1"), alpha = 0.50, cat.col = rep('black', 2), 
             col = 'black', cex = 1.5, fontfamily = 'serif', 
             cat.cex = 1.5, cat.fontfamily = 'serif',filename = NULL)

pdf('KEGG_mb_up_venn.pdf', width = 10, height = 10)
grid.draw(p_venn_up)
dev.off()



venn_list_down <- list(metagenomic = KEGG_result_down$Description, metabolite = result_all_down_pathway)

p_venn_down <- venn.diagram(venn_list_down, imagetype = 'png', 
                          fill = c("dodgerblue", "goldenrod1"), alpha = 0.50, cat.col = rep('black', 2), 
                          col = 'black', cex = 1.5, fontfamily = 'serif', 
                          cat.cex = 1.5, cat.fontfamily = 'serif',filename = NULL)

pdf('KEGG_mb_down_venn.pdf', width = 10, height = 10)
grid.draw(p_venn_down)
dev.off()


#继续以上述4个分组为例，组间交集元素获得
inter <- get.venn.partitions(venn_list_up)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(3, 4)], 'up_venn_inter.txt', row.names = FALSE, sep = '\t', quote = FALSE)


inter <- get.venn.partitions(venn_list_down)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(3, 4)], 'down_venn_inter.txt', row.names = FALSE, sep = '\t', quote = FALSE)








