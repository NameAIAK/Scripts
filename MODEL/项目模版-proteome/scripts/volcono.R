rm(list=ls())
####load packages#####
library(EnhancedVolcano)
library(limma)
library(edgeR)
library(readxl)
library(dplyr)
library(metPath)
library(ggplot2)
library(org.Hs.eg.db)
library(org.Mm.eg.db) 
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggvenn)
library(VennDiagram)

P.Value <- 0.05
logFC <- 0


####M-C#####
norm_counts <- read_excel('Yang-20241207-Protein.xlsx')

ACC_EN_GENE <- norm_counts[,c('Accession','Ensembl Gene ID','Gene Symbol')]

norm_counts1<- norm_counts[,c('Accession',
                              'Abundances (Normalized): F1: Sample, 1',
                              'Abundances (Normalized): F2: Sample, 1',
                              'Abundances (Normalized): F3: Sample, 1',
                              'Abundances (Normalized): F4: Sample, 2',
                              'Abundances (Normalized): F5: Sample, 2',
                              'Abundances (Normalized): F6: Sample, 2'
)]
accession <- norm_counts1$Accession
norm_counts1 <- norm_counts1[,-1]
norm_counts1[is.na(norm_counts1)] <- 0
rownames(norm_counts1) <- accession

dge <- DGEList(counts=norm_counts1)
group.list=c(rep('C',3),rep('M',3))
group.list=factor(group.list)
group.list
group.list=relevel(group.list,ref = "C")

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
# DEG_MC=topTable(fit2, coef=1, adjust="BH")  
DEG_MC = na.omit(DEG)   
write.csv(DEG_MC,file = 'proteome_MC.csv')
colnames(DEG_MC)
DEG_MC$change = ifelse(DEG_MC$P.Value< P.Value & abs(DEG_MC$logFC) > logFC, 
                   ifelse(DEG_MC$logFC> logFC ,'Up','Down'),
                   'Stable')
p <- ggplot(
  # 数据、映射、颜色
  DEG_MC, aes(x = logFC, y = -log10(P.Value), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="logFC",
       y="-log10(P.Value)")+
  scale_x_continuous(limits = c(-5, 5))+
  # scale_y_continuous(limits = c(0,30))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('M v.s. C')
p
# 将火山图另存为SVG文件
ggsave(
  plot = p,
  "MC.p_volcano.pdf",
  height = 6,
  width = 6
)

df_up <- DEG_MC[DEG_MC$logFC > logFC & DEG_MC[, 'P.Value'] < P.Value,]
df_up <- df_up[complete.cases(df_up), ] 
df_upf <- rownames(df_up)
# df_upf

df_up <- df_up %>%  
  mutate(change = "UP")  
write.table(df_up, file = "MC.up.csv", sep = ",", row.names = TRUE, quote = FALSE, na = "")

df_down <- DEG_MC[DEG_MC$logFC < -logFC & DEG_MC[, 'P.Value'] < P.Value,]
df_down <- df_down[complete.cases(df_down), ] 
df_downf <- rownames(df_down)
# df_downf
df_down <- df_down %>%  
  mutate(change = "DOWN")
write.table(df_down, file = "MC.down.csv", sep = ",", row.names = TRUE, quote = FALSE, na = "")

change <- bind_rows(df_up,df_down)
write.table(change, file = "MC.change.csv", sep = ",", row.names = TRUE, quote = FALSE, na = "")


up_acc <- as.data.frame(df_upf)
colnames(up_acc) <- 'Accession'
up_acc <- merge(up_acc,ACC_EN_GENE,by='Accession')
up_acc <- up_acc[!(is.na(up_acc$`Ensembl Gene ID`) & is.na(up_acc$`Gene Symbol`)),]
MC_UP <- up_acc
write.csv(MC_UP,file = 'P_MC_UP.csv',row.names = FALSE)
# table(is.na(up_acc$`Ensembl Gene ID`))
# table(is.na(up_acc$`Gene Symbol`))
down_acc <- as.data.frame(df_downf)
colnames(down_acc) <- 'Accession'
down_acc <- merge(down_acc,ACC_EN_GENE,by='Accession')
down_acc <- down_acc[!(is.na(down_acc$`Ensembl Gene ID`) & is.na(down_acc$`Gene Symbol`)),]
MC_DOWN <- down_acc
write.csv(MC_DOWN,file = 'P_MC_DOWN.csv',row.names = FALSE)
####MC-enrichment####
####MC-UP####
gs = bitr(unique(up_acc$`Gene Symbol`), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
####MC-GO_up####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_up <- data.frame(ego.bp@result)
ego.bp_up <- ego.bp_up[ego.bp_up$pvalue<0.05,]
# ego.bp_up <- ego.bp_up[c(1:10),]
# # ego.bp_up$GeneRatio_numeric <- sapply(ego.bp_up$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_up[order(ego.bp_up$GeneRatio_numeric),]$Description

ego.bp_up$Description<- factor(ego.bp_up$Description,levels =ego.bp_up[order(ego.bp_up$Count),]$Description)
write.csv(file="MC_ego.bp_up_p0.05.csv",data.frame(ego.bp_up),row.names=F)

pdf(file="MC_ego.bp_up_bar.pdf",width = 6,height = 7)
ggplot(ego.bp_up[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="MC GO Enrichment UP")+
  theme_bw()
dev.off()


####MC-kegg_up####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'mmu',pvalueCutoff = 0.05)
kk_up <- data.frame(kk@result)
kk_up <- kk_up[kk_up$pvalue<0.05,]
kk_up <- kk_up %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))
# kk_up <- ego.bp_up[c(1:10),]
# # ego.bp_up$GeneRatio_numeric <- sapply(ego.bp_up$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_up[order(ego.bp_up$GeneRatio_numeric),]$Description

kk_up$Description<- factor(kk_up$Description,levels =kk_up[order(kk_up$Count),]$Description)
write.csv(file="MC_kk_up_p0.05.csv",data.frame(kk_up),row.names=F)

pdf(file="MC_kk_up_bar.pdf",width = 7,height = 7)

ggplot(data = kk_up[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "MC_KEGG_UP") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()

####MC-DOWN####
gs = bitr(unique(down_acc$`Gene Symbol`), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
####MC-GO_down####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_down <- data.frame(ego.bp@result)
ego.bp_down <- ego.bp_down[ego.bp_down$pvalue<0.05,]
# ego.bp_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description

ego.bp_down$Description<- factor(ego.bp_down$Description,levels =ego.bp_down[order(ego.bp_down$Count),]$Description)
write.csv(file="MC_ego.bp_down_p0.05.csv",data.frame(ego.bp_down),row.names=F)

pdf(file="MC_ego.bp_down_bar.pdf",width = 8.5,height = 7)
ggplot(ego.bp_down[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="MC GO Enrichment DOWN")+
  theme_bw()
dev.off()


####MC-kegg_down####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'mmu',pvalueCutoff = 0.05)
kk_down <- data.frame(kk@result)
kk_down <- kk_down[kk_down$pvalue<0.05,]
kk_down <- kk_down %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))
# kk_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description

kk_down$Description<- factor(kk_down$Description,levels =kk_down[order(kk_down$Count),]$Description)
write.csv(file="MC_kk_down_p0.05.csv",data.frame(kk_down),row.names=F)


pdf(file="MC_kk_down_bar.pdf",width = 7,height = 7)

ggplot(data = kk_down[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "MC_KEGG_down") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()



####T-M####
norm_counts2 <- norm_counts[,c('Accession',
                              'Abundances (Normalized): F4: Sample, 2',
                              'Abundances (Normalized): F5: Sample, 2',
                              'Abundances (Normalized): F6: Sample, 2',
                              'Abundances (Normalized): F7: Sample, 4',
                              'Abundances (Normalized): F8: Sample, 4',
                              'Abundances (Normalized): F9: Sample, 4'
)]
accession <- norm_counts2$Accession
norm_counts2 <- norm_counts2[,-1]
norm_counts2[is.na(norm_counts2)] <- 0
rownames(norm_counts2) <- accession

dge <- DGEList(counts=norm_counts2)
group.list=c(rep('M',3),rep('T',3))
group.list=factor(group.list)
group.list
group.list=relevel(group.list,ref = "M")

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

# DEG2=topTable(fit2, coef=1, adjust="BH")  
DEG_TM = na.omit(DEG)   
write.csv(DEG_TM,file = 'proteome_TM.csv')

DEG_TM$change = ifelse(DEG_TM$P.Value< P.Value & abs(DEG_TM$logFC) > logFC, 
                       ifelse(DEG_TM$logFC> logFC ,'Up','Down'),
                       'Stable')
p <- ggplot(
  # 数据、映射、颜色
  DEG_TM, aes(x = logFC, y = -log10(P.Value), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept =-log10(P.Value),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="logFC",
       y="-log10(P.Value)")+
  scale_x_continuous(limits = c(-5, 5))+
  # scale_y_continuous(limits = c(0,30))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('T v.s. M')
p
# 将火山图另存为SVG文件
ggsave(
  plot = p,
  "TM.p_volcano.pdf",
  height = 6,
  width = 6
)

df_up <- DEG_TM[DEG_TM$logFC > logFC & DEG_TM[, 'P.Value'] < P.Value,]
df_up <- df_up[complete.cases(df_up), ] 
df_upf <- rownames(df_up)
# df_upf

df_up <- df_up %>%  
  mutate(change = "UP")  
write.table(df_up, file = "TM.up.csv", sep = ",", row.names = TRUE, quote = FALSE, na = "")

df_down <- DEG_TM[DEG_TM$logFC < -logFC & DEG_TM[, 'P.Value'] < P.Value,]
df_down <- df_down[complete.cases(df_down), ] 
df_downf <- rownames(df_down)
# df_downf
df_down <- df_down %>%  
  mutate(change = "DOWN")
write.table(df_down, file = "TM.down.csv", sep = ",", row.names = TRUE, quote = FALSE, na = "")

change <- bind_rows(df_up,df_down)
write.table(change, file = "TM.change.csv", sep = ",", row.names = TRUE, quote = FALSE, na = "")

up_acc <- as.data.frame(df_upf)
colnames(up_acc) <- 'Accession'
up_acc <- merge(up_acc,ACC_EN_GENE,by='Accession')

up_acc <- up_acc[!(is.na(up_acc$`Ensembl Gene ID`) & is.na(up_acc$`Gene Symbol`)),]
TM_UP <- up_acc
write.csv(TM_UP,file = 'P_TM_UP.csv',row.names = FALSE)
# table(is.na(up_acc$`Ensembl Gene ID`))
# table(is.na(up_acc$`Gene Symbol`))

down_acc <- as.data.frame(df_downf)
colnames(down_acc) <- 'Accession'
down_acc <- merge(down_acc,ACC_EN_GENE,by='Accession')

down_acc <- down_acc[!(is.na(down_acc$`Ensembl Gene ID`) & is.na(down_acc$`Gene Symbol`)),]
TM_DOWN <- down_acc
write.csv(TM_DOWN,file = 'P_TM_DOWN.csv',row.names = FALSE)
#######TM-enrichment#######
####TM-UP####
gs = bitr(unique(up_acc$`Gene Symbol`), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
####TM-GO_up####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_up <- data.frame(ego.bp@result)
ego.bp_up <- ego.bp_up[ego.bp_up$pvalue<0.05,]
# ego.bp_up <- ego.bp_up[c(1:10),]
# # ego.bp_up$GeneRatio_numeric <- sapply(ego.bp_up$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_up[order(ego.bp_up$GeneRatio_numeric),]$Description

ego.bp_up$Description<- factor(ego.bp_up$Description,levels =ego.bp_up[order(ego.bp_up$Count),]$Description)
write.csv(file="TM_ego.bp_up_p0.05.csv",data.frame(ego.bp_up),row.names=F)

pdf(file="TM_ego.bp_up_bar.pdf",width = 6,height = 7)
ggplot(ego.bp_up[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="TM GO Enrichment UP")+
  theme_bw()
dev.off()


####TM-kegg_up####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'mmu',pvalueCutoff = 0.05)
kk_up <- data.frame(kk@result)
kk_up <- kk_up[kk_up$pvalue<0.05,]
kk_up <- kk_up %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))

# kk_up <- ego.bp_up[c(1:10),]
# # ego.bp_up$GeneRatio_numeric <- sapply(ego.bp_up$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_up[order(ego.bp_up$GeneRatio_numeric),]$Description

kk_up$Description<- factor(kk_up$Description,levels =kk_up[order(kk_up$Count),]$Description)
write.csv(file="TM_kk_up_p0.05.csv",data.frame(kk_up),row.names=F)


pdf(file="TM_kk_up_bar.pdf",width = 7,height = 7)

ggplot(data = kk_up[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "TM_KEGG_UP") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()

####TM-DOWN####
gs = bitr(unique(down_acc$`Gene Symbol`), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
####TM-GO_down#####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_down <- data.frame(ego.bp@result)
ego.bp_down <- ego.bp_down[ego.bp_down$pvalue<0.05,]
# ego.bp_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description

ego.bp_down$Description<- factor(ego.bp_down$Description,levels =ego.bp_down[order(ego.bp_down$Count),]$Description)
write.csv(file="TM_ego.bp_down_p0.05.csv",data.frame(ego.bp_down),row.names=F)

pdf(file="TM_ego.bp_down_bar.pdf",width = 6,height = 7)
ggplot(ego.bp_down[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="TM GO Enrichment DOWN")+
  theme_bw()
dev.off()


####TM-kegg_down####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'mmu',pvalueCutoff = 0.05)
kk_down <- data.frame(kk@result)
kk_down <- kk_down[kk_down$pvalue<0.05,]
kk_down <- kk_down %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))
# kk_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description

kk_down$Description<- factor(kk_down$Description,levels =kk_down[order(kk_down$Count),]$Description)
write.csv(file="TM_kk_down_p0.05.csv",data.frame(kk_down),row.names=F)


pdf(file="TM_kk_down_bar.pdf",width = 7,height = 7)

ggplot(data = kk_down[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "TM_KEGG_down") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()

####MC/TM inner gene####
MCUPinnerTMDOWN<- list('M vs C'=MC_UP$`Gene Symbol`,'T vs M'=TM_DOWN$`Gene Symbol`)
MCDOWNinnerTMUP<-list('M vs C'=MC_DOWN$`Gene Symbol`,'T vs M'=TM_UP$`Gene Symbol`)

####MCUPinnerTMDOWN####
p_MCUPinnerTMDOWN <- ggvenn(MCUPinnerTMDOWN,c('M vs C', 'T vs M'),show_percentage = T,
       stroke_color = "white",
       fill_color = c("#BC3C29FF","#0072B5FF","#E18727FF","#FFDC91FF",
                      "#EE4C97FF","#20854EFF","#7876B1FF","#6F99ADFF" ),
       set_name_color =c("#E41A1C","#1E90FF"))
ggsave(p_MCUPinnerTMDOWN,filename = 'MCUPinnerTMDOWN.venn.pdf',width = 6,height = 4)

inner1 <- get.venn.partitions(MCUPinnerTMDOWN)
inner1_gene <- data.frame(inner1[[1,'..values..']])
colnames(inner1_gene) <- 'MCUPinnerTMDOWN_gene'
write.table(inner1_gene, 'MCUPinnerTMDOWN_gene_inter.csv', row.names = FALSE, col.names = TRUE,sep = ',', quote = FALSE)

gs = bitr(inner1_gene$MCUPinnerTMDOWN_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
####inner1GO#####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_inner1 <- data.frame(ego.bp@result)
ego.bp_inner1 <- ego.bp_inner1[ego.bp_inner1$pvalue<0.05,]
# ego.bp_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description
ego.bp_inner1$Description<- factor(ego.bp_inner1$Description,levels =ego.bp_inner1[order(ego.bp_inner1$Count),]$Description)
write.csv(file="MCUPinnerTMDOWN_GO_p0.05.csv",data.frame(ego.bp_inner1),row.names=F)

pdf(file="MCUPinnerTMDOWN_GO.pdf",width = 7.5,height = 7)
ggplot(ego.bp_inner1[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="MCUPinnerTMDOWN_GO")+
  theme_bw()
dev.off()

####inner1KEGG####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'mmu',pvalueCutoff = 0.05)
kk_inner1 <- data.frame(kk@result)
kk_inner1 <- kk_inner1[kk_inner1$pvalue<0.05,]
kk_inner1 <- kk_inner1 %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))
# kk_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description

kk_inner1$Description<- factor(kk_inner1$Description,levels =kk_inner1[order(kk_inner1$Count),]$Description)
write.csv(file="MCUPinnerTMDOWN_KEGG_p0.05.csv",data.frame(kk_inner1),row.names=F)


pdf(file="MCUPinnerTMDOWN_KEGG_bar.pdf",width = 7,height = 7)

ggplot(data = kk_inner1[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "MCUPinnerTMDOWN_KEGG") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()

####MCDOWNinnerTMUP####
p_MCDOWNinnerTMUP <- ggvenn(MCDOWNinnerTMUP,c('M vs C', 'T vs M'),show_percentage = T,
                            stroke_color = "white",
                            fill_color = c("#BC3C29FF","#0072B5FF","#E18727FF","#FFDC91FF",
                                           "#EE4C97FF","#20854EFF","#7876B1FF","#6F99ADFF" ),
                            set_name_color =c("#E41A1C","#1E90FF"))
ggsave(p_MCDOWNinnerTMUP,filename = 'MCDOWNinnerTMUP.venn.pdf',width = 6,height = 4)

inner2 <- get.venn.partitions(MCDOWNinnerTMUP)
inner2_gene <- data.frame(inner2[[1,'..values..']])
colnames(inner2_gene) <- 'MCDOWNinnerTMUP_gene'
write.table(inner2_gene, 'MCDOWNinnerTMUP_gene_inter.csv', row.names = FALSE, col.names = TRUE,sep = ',', quote = FALSE)

gs = bitr(inner2_gene$MCDOWNinnerTMUP_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
####inner2GO#####
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_inner2 <- data.frame(ego.bp@result)
ego.bp_inner2 <- ego.bp_inner2[ego.bp_inner2$pvalue<0.05,]
# ego.bp_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description
ego.bp_inner2$Description<- factor(ego.bp_inner2$Description,levels =ego.bp_inner2[order(ego.bp_inner2$Count),]$Description)
write.csv(file="MCDOWNinnerTMUP_GO_p0.05.csv",data.frame(ego.bp_inner2),row.names=F)

pdf(file="MCDOWNinnerTMUP_GO.pdf",width =6.5,height = 7)
ggplot(ego.bp_inner2[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="MCDOWNinnerTMUP_GO")+
  theme_bw()
dev.off()

####inner2KEGG####
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'mmu',pvalueCutoff = 0.05)
kk_inner2 <- data.frame(kk@result)
kk_inner2 <- kk_inner2[kk_inner2$pvalue<0.05,]
kk_inner2 <- kk_inner2 %>%
  filter(!grepl("disease", Description, ignore.case = TRUE))
# kk_down <- ego.bp_down[c(1:10),]
# # ego.bp_down$GeneRatio_numeric <- sapply(ego.bp_down$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_down[order(ego.bp_down$GeneRatio_numeric),]$Description

kk_inner2$Description<- factor(kk_inner2$Description,levels =kk_inner2[order(kk_inner2$Count),]$Description)
write.csv(file="MCDOWNinnerTMUP_KEGG_p0.05.csv",data.frame(kk_inner2),row.names=F)

pdf(file="MCDOWNinnerTMUP_KEGG_bar.pdf",width = 7,height = 7)
ggplot(data = kk_inner2[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "MCDOWNinnerTMUP_KEGG") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()
