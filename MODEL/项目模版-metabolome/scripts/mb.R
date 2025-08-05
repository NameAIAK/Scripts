
library(clusterProfiler)
library(readxl)
library(data.table)
library("DESeq2")
library("EnhancedVolcano")
library('dplyr') 


metabolite_row <- read_excel("metabolite_data.xlsx")
ncol_metabolite <- ncol(metabolite_row)
metabolite <- metabolite_row[,c(2,11,17:ncol_metabolite)]

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
write.csv(DEG,file = 'mb.P_C v.s. O.csv')

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
write.csv(DEG,file = 'mb N_C v.s. O.csv')

neg_diff_mb <- DEG[DEG[, "P.Value"] < 0.05,]
variables_neg <- rownames(neg_diff_mb)

n_mb_up <- neg_diff_mb[neg_diff_mb$logFC>0,]
variables_neg_up<- rownames(n_mb_up)

n_mb_down <- neg_diff_mb[neg_diff_mb$logFC<0,]
variables_neg_down <- rownames(n_mb_down)


df_neg <- df_mb_neg[,c(variables_neg)]


keggannotation <- read.csv('keggannotation.csv')
keggannotation <- keggannotation[,c(2,3)]
# write.csv(keggannotation,file='keggannotation.csv')

allkeggid <- metabolite$KEGG
allkeggid <- data.frame(allkeggid)
colnames(allkeggid) <- c('ID')
total <- right_join(keggannotation,allkeggid,by="ID") %>% select(2,1)


mb_result_down <- clusterProfiler::enricher(gene = variables_pos_down,TERM2GENE = total,minGSSize = 1,pvalueCutoff = 1,qvalueCutoff = 1)
mb_result_down <- mb_result_down@result
mb_result_down <- mb_result_down[mb_result_down$pvalue<0.05,]




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


