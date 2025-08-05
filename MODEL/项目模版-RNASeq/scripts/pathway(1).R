
######使用clusterProfiler包，GO和KEGG富集分析
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db) ##加载人类
keytypes(org.Mm.eg.db)
DEG<- c(RHM[RHM$change=='Down',]$gene)

gs = bitr(DEG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(gs)
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
write.csv(ego.bp,file = 'D:\\data\\洪祝平\\new结果(-A)\\GO.csv',row.names = F)
p1<- print(dotplot(ego.bp, showCategory=10,title="GO"))

kk <- enrichKEGG(gene= gs$ENTREZID, organism     = 'mmu',pvalueCutoff = 0.05)
write.csv(kk,file = 'KEGG_upno.csv',row.names = F)
p2<- print(dotplot(kk, showCategory=10,title="KEGG_upno")) #展示前十个条???
p1 | p2

colnames(dat)
DEG<- dat[intersect(which(abs(dat$logFC.YJSTW16h)>1),which(dat$FDR.y.1<0.05)),1]


pdf('./DEG功能富集/YJSTW16h.pdf',12,8)
p1 | p2
dev.off()

###########GSEA
library(msigdbr)
library(fgsea)
library(dplyr)
library(tibble)

MC<- read.table('RHvsM_2.txt')

head(MC)
MC$gene<- rownames(MC)
cluster.genes<- MC %>% arrange(desc(logFC)) %>% dplyr::select(gene,logFC) #基因按logFC排序
ranks<- deframe(cluster.genes)

genelist<- list.files('./胆汁酸')

#Mm_m5 <- msigdbr(species = "Mus musculus", category = "m5")## 定义基因集，选取C2
#fgsea_sets = mdb_c2 [grep("^GO",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets<- list()
for (i in 1:length(genelist)) {
  dir<- paste0('./胆汁酸/',genelist[i])
  temp<- read.gmt(dir)
  fgsea_sets[[names(table(temp$term))]]<- temp$gene
  
}

fgseaRes<- fgsea::fgsea(fgsea_sets, stats = ranks) #运行fgsea
library(fgsea)


plotEnrichment(fgsea_sets[[2]],ranks) + 
  labs(title=names(fgsea_sets)[2]) #对某???特定通路分析


####2.clusterProfiler包分???
#报错
library(stringr)
library(enrichplot)
library(ggplot2)
library(org.Mm.eg.db)
geneset <- read.gmt("./胆汁酸/REACTOME_SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS.v2024.1.Mm.gmt")  
geneset<- geneset[geneset$term %in% c(names(fgsea_sets)),]
head(geneset)

gs <-bitr(MC$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
markers1<-cbind(MC[gs[,1],],gs)
geneList = markers1$logFC
names(geneList) = markers1$ENTREZID
geneList = sort(geneList,decreasing = T)

gs <-bitr(geneset$gene, fromType="SYMBOL", toType="ENTREZID", 
          OrgDb="org.Mm.eg.db")
geneset<-geneset[geneset$gene %in% gs$SYMBOL,]
geneset$gene<- gs$ENTREZID

egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F,
             pAdjustMethod = 'none',pvalueCutoff = 1)
egmt1<- setReadable(egmt,OrgDb=org.Mm.eg.db, keyType = "ENTREZID")
y=data.frame(egmt1)


library(enrichplot)

gseaplot2(egmt, geneSetID =1, pvalue_table = TRUE,
          title = egmt@result$ID[1])

egmt@result$NES

# 感兴趣的 Gene
symbol <- egmt[[egmt$ID[[1]]]]
head(symbol)

g<- egmt1[[egmt1$ID[[1]]]]
p+geom_gsea_gene(g)

features <- g

data<- read.csv('CMM.csv')
rownames(data)<- data[,1]
data<- data[,-1]
dd<- data[c(up$`M vs C`,up$`M vs A`,up$`M vs B`,
            down$`M vs C`,down$`M vs A`,down$`M vs B`),]
pheatmap::pheatmap(dd,scale = 'row')


#########GSVA分析
library(ggplot2)
library(dplyr)
library(msigdbr)
library(GSVA)
library(pheatmap)
library(patchwork)
genelist<- list.files('./GENESETS')

#Mm_m5 <- msigdbr(species = "Mus musculus", category = "m5")## 定义基因集，选取C2
#fgsea_sets = mdb_c2 [grep("^GO",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets<- list()
for (i in 1:length(genelist)) {
  dir<- paste0('./GENESETS/',genelist[i])
  temp<- read.gmt(dir)
  fgsea_sets[[names(table(temp$term))]]<- temp$gene
  
}


expr=as.matrix(norm_counts) 
kegg <- gsva(expr, fgsea_sets, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
pheatmap(kegg)#绘制热图


library(pheatmap)
p <- pheatmap(kegg[c(3,4,10,11,12,13,14),],
              cluster_rows = F,
              cluster_cols = F,
              show_rownames = T,
              show_colnames = T,
              color =colorRampPalette(c("blue", "white","red"))(100))
print(p)

######################
library(forcats)
library(dittoSeq)

kegg1<- data.frame(BioPlanet_2019_table)
head(kegg1)
kegg1<- kegg1[1:10,]
kegg1$Overlap
kegg1$Overlap<- c(4/24,3/17,4/47,2/5,2/6,10/432,2/8,3/32,12/603,3/34)
p <- ggplot(kegg1,aes(Overlap,Term))+
  geom_bar(stat = "identity", aes(fill=P.value))+

  labs(y='Pathway',x='Gene Ratio')+
  theme_classic()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text = element_text(colour = 'black', size = 10),
        plot.margin = margin(0,0,0,-0.05, "cm"))+
  scale_fill_manual(values = dittoColors())
p

ggplot(kegg1,aes(x=Overlap,
                 y=reorder(Term,Overlap)))+
  geom_col(aes(fill=P.value))+
  scale_fill_distiller(palette = "RdPu",direction = 1) +

  labs(x="Gene Ratio",y='Pathway')+
  theme_bw()+
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,0.65))+
  theme(panel.grid = element_blank())
################
kegg1<- read.csv('GO_downup.csv')
head(kegg1)
which(kegg1$Description=='bile acid biosynthetic process')
#4,18,37,49,127
#1,2,76,99,118

#3,6,11,14,19
#1,2,7,12,16

top10<- kegg1[c(4,18,37,49,127),]
top10$GeneRatio
top10$GeneRatio<- c(41/346,23/346,8/346,7/346,4/346)

ggplot(top10,aes(x=GeneRatio,
                 y=reorder(Description,GeneRatio)))+
  geom_col(aes(fill=p.adjust))+
  scale_fill_distiller(palette = "YlOrRd",direction = 1) +
  
  labs(x="Gene Ratio",y='Pathway')+
  theme_bw()+
  theme(panel.grid = element_blank())

