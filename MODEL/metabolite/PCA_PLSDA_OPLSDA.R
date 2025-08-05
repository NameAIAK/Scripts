# install ropls
if (F) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("ropls")
}
 
# load  packages
library(ropls)
library(ggplot2)
library(ggsci)
library(Cairo)
library(tidyverse)
library(extrafont)
loadfonts()

# load data
data(sacurine)
names(sacurine)
 
# view data information
attach(sacurine)
strF(dataMatrix)
strF(sampleMetadata)
strF(variableMetadata)

# PCA analysis
pca = opls(dataMatrix)
genderFc = sampleMetadata[, "gender"]
 
pdf(file = 'figures/PCA.pdf', width = 5, height = 5)
plot(pca, typeVc = "x-score",
     parAsColFcVn = genderFc, parEllipsesL = TRUE)
dev.off()

# PLSDA analysis
plsda = opls(dataMatrix,genderFc)
 
# sample scores plot
sample.score = plsda@scoreMN %>% 
  as.data.frame() %>%
  mutate(gender = sacurine[["sampleMetadata"]][["gender"]])
  
p1 = ggplot(sample.score, aes(p1, p2, color = gender)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_point() +
  geom_point(aes(-10,-10), color = 'white') +
  labs(x = 'P1(10.0%)',y = 'P2(9%)') +
  stat_ellipse(level = 0.95, linetype = 'solid', 
               size = 1, show.legend = FALSE) +
  scale_color_manual(values = c('#008000','#FFA74F')) +
  theme_bw() +
  theme(legend.position = c(0.9,0.8),
    legend.text = element_text(color = 'black',size = 12, family = 'Arial', face = 'plain'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
    axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
    axis.ticks = element_line(color = 'black'))
ggsave(p1, filename = 'figures/pls.pdf', 
       width = 5, height = 5, device = cairo_pdf)



# VIP scores plot
vip.score = as.data.frame(plsda@vipVn)
colnames(vip.score) = 'vip'
vip.score$metabolites = rownames(vip.score)
vip.score = vip.score[order(-vip.score$vip),]
vip.score$metabolites = factor(vip.score$metabolites,
                               levels = vip.score$metabolites)
 
loading.score = plsda@loadingMN %>% as.data.frame()
loading.score$metabolites = rownames(loading.score)
 
all.score = merge(vip.score, loading.score, by = 'metabolites')
 
all.score$cat = paste('A',1:nrow(all.score), sep = '')
 
p2 = ggplot(all.score[all.score$vip >= 1,], aes(cat, vip)) +
  geom_segment(aes(x = cat, xend = cat,
                   y = 0, yend = vip)) +
  geom_point(shape = 21, size = 5, color = '#008000' ,fill = '#008000') +
  geom_point(aes(1,2.5), color = 'white') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'VIP value') +
  theme_bw() +
  theme(legend.position = 'none',
        legend.text = element_text(color = 'black',size = 12, family = 'Arial', face = 'plain'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.ticks = element_line(color = 'black'),
        axis.ticks.x = element_blank())
ggsave(p2, filename = 'figures/pls_VIP.pdf', 
       width = 8, height = 5, device = cairo_pdf)

# OPLS-DA analysis
oplsda = opls(dataMatrix, genderFc, predI = 1, orthoI = NA)
 
# sample scores plot
sample.score = oplsda@scoreMN %>% 
  as.data.frame() %>%
  mutate(gender = sacurine[["sampleMetadata"]][["gender"]],
         o1 = oplsda@orthoScoreMN[,1])
 
p3 = ggplot(sample.score, aes(p1, o1, color = gender)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_point() +
  #geom_point(aes(-10,-10), color = 'white') +
  labs(x = 'P1(5.0%)',y = 'to1') +
  stat_ellipse(level = 0.95, linetype = 'solid', 
               size = 1, show.legend = FALSE) +
  scale_color_manual(values = c('#008000','#FFA74F')) +
  theme_bw() +
  theme(legend.position = c(0.1,0.85),
        legend.title = element_blank(),
        legend.text = element_text(color = 'black',size = 12, family = 'Arial', face = 'plain'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.ticks = element_line(color = 'black'))
ggsave(p3, filename = 'figures/opls.pdf', 
       width = 5, height = 5, device = cairo_pdf)

# VIP scores plot
vip.score = as.data.frame(oplsda@vipVn)
colnames(vip.score) = 'vip'
vip.score$metabolites = rownames(vip.score)
vip.score = vip.score[order(-vip.score$vip),]
vip.score$metabolites = factor(vip.score$metabolites,
                               levels = vip.score$metabolites)
 
loading.score = oplsda@loadingMN %>% as.data.frame()
loading.score$metabolites = rownames(loading.score)
 
all.score = merge(vip.score, loading.score, by = 'metabolites')
 
all.score$cat = paste('A',1:nrow(all.score), sep = '')
 
p4 = ggplot(all.score[all.score$vip >= 1,], aes(cat, vip)) +
  geom_segment(aes(x = cat, xend = cat,
                   y = 0, yend = vip)) +
  geom_point(shape = 21, size = 5, color = '#008000' ,fill = '#008000') +
  geom_point(aes(1,2.5), color = 'white') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'VIP value') +
  theme_bw() +
  theme(legend.position = 'none',
        legend.text = element_text(color = 'black',size = 12, family = 'Arial', face = 'plain'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.ticks = element_line(color = 'black'),
        axis.ticks.x = element_blank())
p4
ggsave(p4, filename = 'figures/opls_VIP.pdf', 
       width = 8, height = 5, device = cairo_pdf)



# model training
oplsda.2 = opls(dataMatrix, genderFc, predI = 1, orthoI = NA,subset = "odd") 
OPLS-DA
# 模型在训练集上的准确率
trainVi = getSubsetVi(oplsda.2)
tab = table(genderFc[trainVi], fitted(oplsda.2))
print(paste('模型准确率：',round(sum(diag(tab))/sum(tab)*100, 2),'%', sep = ''))
# model on test data
tab2 = table(genderFc[-trainVi],predict(oplsda.2, dataMatrix[-trainVi, ]))
print(paste('模型准确率：',round(sum(diag(tab2))/sum(tab2)*100, 2),'%', sep = ''))


