### 单细胞分析相关的包
library(AUCell)
library(CellChat)
library(GSVA)
library(RColorBrewer)
library(Seurat)
library(SeuratData)
library(clusterProfiler)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(harmony)
library(limma)
library(monocle)
library(patchwork)
library(pheatmap)
library(progeny)
library(readr)
library(tibble)
library(tidyr)
library(tidyverse)
# monocle3不能与monocle（monocle2）同时加载
library(monocle3)

### 一些单细胞数据处理的小方法
# 查看对象的基因数，细胞数
dim(scRNA)  # 返回基因（行）和细胞（列）的数量
head(scRNA)  # 查看前6行的数据
scRNA[1:5, 1:5]  # 查看前5行和前5列的数据
str(scRNA)      # 查看数据的结构信息

Idents(scRNA)  # 查看所有细胞的聚类标签，基本不使用，打印所有细胞的标签，屏幕打印过多
table(Idents(scRNA))  # 统计每个聚类的细胞数量
head(Idents(obj))   # 显示前几个细胞的分类信息。
Idents(obj)='orig.ident'   #指定使用的分类信息
scRNA@meta.database #存储样本信息，分类信息
table(scRNA@meta.database$orig.ident)  #对分类信息进行统计 

active.assay#查看当前使用的assays
active.ident#查看当前的使用分群方式(可使用 levels 函数)

DefaultAssay(scRNA)   # 查看当前激活的 assay（默认是 "RNA"）
Assays(obj) # 查看所有可用的 assay
DefaultAssay(scRNA) <- "RNA"    #指定使用的数据
DefaultDimReduc(obj)    # 查看当前激活的降维结果（默认是最后一个运行的降维方法）
Reductions(obj) # 查看所有可用的降维结果

saveRDS(scRNA, file = "scRNA_processed.rds")  # 保存处理后的数据
scRNA <- readRDS("scRNA_processed.rds")       # 加载保存的数据

##使用频率低
cell_counts <- colSums(scRNA)  # 计算每个细胞的总表达量（UMI计数）
head(cell_counts)
rownames(scRNA)  # 查看所有基因的名称
colnames(scRNA)  # 查看所有细胞的名称
summary(scRNA)  # 查看数据的基本统计信息
# 提取表达矩阵
expression_matrix <- GetAssayData(obj, slot = "counts")
# 计算每个基因的表达频率（在多少个细胞中表达）
gene_counts <- rowSums(expression_matrix > 0)
head(gene_counts)



rownames(object)#获取全部基因ID
Cells(object)#获取整个object的细胞ID
colnames(object)#获取整个object的细胞ID
Whichcells(object,idents =c(1,2))#按照idents获取部分细胞ID
Whichcells(object,expression=gene1 >1)#按照基因表达获取部分细胞ID
Whichcells(object,expression =genel >1,slot = "counts")#按照基因表达获取部分细胞ID


# 提取包含部分细胞的对象
cells=Whichcells(object,idents =1)#提取细胞ID
subset(x=object,cells = cells)#按照细胞ID提取
subset(x=object,idents=c(1,2))#按照idents提取
subset(x=object,idents="cluster")#对细胞簇重新命名后为cluster

subset(object,idents =c(1,2),invert=TRUE)# 想要排除1、2细胞类型
subset(x=object,stim =="CONTROL")#按照meta.data中设置过的stim分组信息提取

subset(x=object, RNA_snn_res.2 == 2)#按照某一个resolution下的分群提取

subset(x=object,gene1 >1)#根据某个基因的表达量来提取
subset(x=object,genel>1,slot ="counts")#根据某个基因的表达量来提取

# 每个聚类细胞数占比
prop.table(table(Idents(object)))
prop.table(table(object$RNA snn res.0.3))

# 计算平均表达量
cluster.averages<-AverageExpression(object)

# 修改聚类后的因子水平
Idents(object)<-factor(Idents(object),levels= c(1,2,3,4,9,8,7,6,5,0))





### 数据合并
## 数据合并method1 分别读取样本的数据-分别创建Seurat对象-合并Seurat object
# 分别读取6个样本的数据   Je  Ji  JT  Se  Si  ST
Je.data <- Read10X(data.dir = "./Je/filtered_feature_bc_matrix") 
Ji.data <- Read10X(data.dir = "./Ji/filtered_feature_bc_matrix") 
JT.data <- Read10X(data.dir = "./JT/filtered_feature_bc_matrix") 
Se.data <- Read10X(data.dir = "./Se/filtered_feature_bc_matrix") 
Si.data <- Read10X(data.dir = "./Si/filtered_feature_bc_matrix") 
ST.data <- Read10X(data.dir = "./ST/filtered_feature_bc_matrix") 
# 分别构建Seurat object
Je.obj <- CreateSeuratObject(counts = Je.data, project = "Je")
Ji.obj <- CreateSeuratObject(counts = Ji.data, project = "Ji")
JT.obj <- CreateSeuratObject(counts = JT.data, project = "JT")
Se.obj <- CreateSeuratObject(counts = Se.data, project = "Se")
Si.obj <- CreateSeuratObject(counts = Si.data, project = "Si")
ST.obj <- CreateSeuratObject(counts = ST.data, project = "ST")
# 合并Seurat object
JS.all <- merge(Je.obj, y = c(Ji.obj, JT.obj,Se.obj,Si.obj,ST.obj), add.cell.ids = c("Je", "Ji", "JT",'Se','Si','ST'), project = "J3S3")
##数据合并method2 
# 导入单样本文件-读取并创建Seurat Object
dir = c(
  "./Je/filtered_feature_bc_matrix",
  "./Ji/filtered_feature_bc_matrix",
  "./JT/filtered_feature_bc_matrix",
  "./Se/filtered_feature_bc_matrix",
  "./Si/filtered_feature_bc_matrix",
  "./ST/filtered_feature_bc_matrix"
)
# 读取并创建Seurat Object
scRNAlist <- list()  # 创建一个空的列表,其中包含三个seurat对象-合并Seurat object
names(dir) = c("Je", "Ji", "JT",'Se','Si','ST')
sample_name = c("Je", "Ji", "JT",'Se','Si','ST')
for (i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = sample_name[i])
}
# 合并Seurat object
scRNA <- merge(scRNAlist[[1]],y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]]))
dim(scRNA)


### 画图

