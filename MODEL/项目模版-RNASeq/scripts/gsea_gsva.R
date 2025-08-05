# https://www.jianshu.com/p/0098baf2df46

1.R包和示例数据
rm(list = ls())
library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)
#devtools::install_github("xjsun1221/tinyarray")
library(tinyarray)
load("exp_Group_deg.Rdata")
ls()
## [1] "deg"   "exp"   "Group"
exp_Group_deg.Rdata里面有三个数据，可以在生信星球公众号后台回复“gsva775”获取,也可以用自己的数据去生成，exp是个芯片表达矩阵，Group是表示分组信息的因子型数据，deg是它的limma差异分析结果。

exp[1:4,1:4]
##        GSM1366348 GSM1366349 GSM1366350 GSM1366351
## RFC2     8.932805   8.679543   8.625015   8.637085
## HSPA6    9.383421   8.605809   9.462774   9.898573
## PAX8     7.916751   8.500635   8.258467   8.553656
## GUCA1A   5.085221   2.414033   1.718570   4.311794
dim(exp)
## [1] 20161    22
Group
##  [1] RA      RA      RA      RA      RA      RA     
##  [7] RA      RA      RA      RA      RA      RA     
## [13] RA      control control control control control
## [19] control control control control
## Levels: control RA
head(deg)
##             logFC   AveExpr         t      P.Value
## HADHA   -1.633518 11.204696 -16.14895 8.628803e-14
## LRRFIP1 -1.713884 11.965832 -15.06767 3.572073e-13
## OTUD4   -2.274611  7.651265 -13.43343 3.618635e-12
## FKBP2   -1.935657  8.628444 -13.04331 6.501205e-12
## ILRUN   -1.893808  8.896532 -12.38173 1.812689e-11
## IGF2R   -1.691681 10.302132 -12.00890 3.290947e-11
##            adj.P.Val        B change
## HADHA   2.358899e-09 20.92048   down
## LRRFIP1 4.882577e-09 19.68466   down
## OTUD4   3.297481e-08 17.61994   down
## FKBP2   3.554534e-08 17.08858   down
## ILRUN   6.607251e-08 16.15073   down
## IGF2R   1.124578e-07 15.60094   down
2.端详msigdbr这个包
1.囊括了11个物种
msigdbr_species()
## # A tibble: 11 x 2
##    species_name             species_common_name      
##    <chr>                    <chr>                    
##  1 Bos taurus               cattle                   
##  2 Caenorhabditis elegans   roundworm                
##  3 Canis lupus familiaris   dog                      
##  4 Danio rerio              zebrafish                
##  5 Drosophila melanogaster  fruit fly                
##  6 Gallus gallus            chicken                  
##  7 Homo sapiens             human                    
##  8 Mus musculus             house mouse              
##  9 Rattus norvegicus        Norway rat               
## 10 Saccharomyces cerevisiae baker's or brewer's yeast
## 11 Sus scrofa               pig
2.看人类有哪些基因集
human <- msigdbr(species = "Homo sapiens")
human[1:4,1:4]
## # A tibble: 4 x 4
##   gs_cat gs_subcat      gs_name        entrez_gene
##   <chr>  <chr>          <chr>                <int>
## 1 C3     MIR:MIR_Legacy AAACCAC_MIR140       10257
## 2 C3     MIR:MIR_Legacy AAACCAC_MIR140       23172
## 3 C3     MIR:MIR_Legacy AAACCAC_MIR140          81
## 4 C3     MIR:MIR_Legacy AAACCAC_MIR140          90
table(human[,1])
## 
##      C1      C2      C3      C4      C5      C6 
##   40056  520069  734948   91173 1223276   30540 
##      C7      C8       H 
##  945462   54466    7321
可以看到，包里和网页上一样，都是这些个collection

H: hallmark gene sets
C1: positional gene sets
C2: curated gene sets
C3: motif gene sets
C4: computational gene sets
C5: GO gene sets
C6: oncogenic signatures
C7: immunologic signatures

3.KEGG和GO
kegg 在C2,GO 在C5。

table(human$gs_subcat)
## 
##                             CGN             CGP 
##         1077845           42544          374969 
##              CM              CP     CP:BIOCARTA 
##           48629            4520            4813 
##         CP:KEGG          CP:PID     CP:REACTOME 
##           12797            8050           88301 
## CP:WIKIPATHWAYS           GO:BP           GO:CC 
##           26619          663648           95709 
##           GO:MF             HPO  MIR:MIR_Legacy 
##          104978          358941           34163 
##       MIR:MIRDB        TFT:GTRD  TFT:TFT_Legacy 
##          372126          173572          155087
KEGG_df = msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_exact_source,gene_symbol)
head(KEGG_df)
## # A tibble: 6 x 2
##   gs_exact_source gene_symbol
##   <chr>           <chr>      
## 1 hsa02010        ABCA1      
## 2 hsa02010        ABCA10     
## 3 hsa02010        ABCA12     
## 4 hsa02010        ABCA13     
## 5 hsa02010        ABCA2      
## 6 hsa02010        ABCA3
# 基因数量
length(unique(KEGG_df$gene_symbol))
## [1] 5245
# 通路数量
length(unique(KEGG_df$gs_exact_source))
## [1] 186
这个数字会随着数据库的更新而发生改变的

GO_df = msigdbr(species = "Homo sapiens",category = "C5") %>% 
  dplyr::select(gene_symbol,gs_exact_source,gs_subcat)
dim(GO_df)
## [1] 1223276       3
GO_df = GO_df[GO_df$gs_subcat!="HPO",]
table(GO_df$gs_subcat)
## 
##  GO:BP  GO:CC  GO:MF 
## 663648  95709 104978
GO_df = GO_df[,c(2,1)]
head(GO_df)
## # A tibble: 6 x 2
##   gs_exact_source gene_symbol
##   <chr>           <chr>      
## 1 GO:0004645      GDPGP1     
## 2 GO:0004645      MTAP       
## 3 GO:0004645      PYGB       
## 4 GO:0004645      PYGL       
## 5 GO:0004645      PYGM       
## 6 GO:0004645      TYMP
# 基因数量
length(unique(GO_df$gene_symbol))
## [1] 19276
# terms数量
length(unique(GO_df$gs_exact_source))
## [1] 10271
可以看到GO包括的基因数量多很多。

3.做GSEA
ge = deg$logFC
names(ge) = rownames(deg)
ge = sort(ge,decreasing = T)
head(ge)
##     OLR1   COL1A1    OLFM4     H3C8   CRISP3  ZFP36L2 
## 3.835314 3.790745 3.427338 3.180043 2.956932 2.909744
TERM2GENE参数可以直接接受上面得到的GO_df作为参数，注意第一列是term第二列是基因哦

head(GO_df)
## # A tibble: 6 x 2
##   gs_exact_source gene_symbol
##   <chr>           <chr>      
## 1 GO:0004645      GDPGP1     
## 2 GO:0004645      MTAP       
## 3 GO:0004645      PYGB       
## 4 GO:0004645      PYGL       
## 5 GO:0004645      PYGM       
## 6 GO:0004645      TYMP
em <- GSEA(ge, TERM2GENE = GO_df)
#画个图来看看
gseaplot2(em, geneSetID = 1, title = em$Description[1])

4.做GSVA
输入数据是表达矩阵,通路和基因间的对应关系需要做成列表的形式，因为gset.idx.list要求是列表。

注意，表达矩阵的行名应与gset.idx.list里的基因名是同一类基因id，比如都是symbol或都是entriz id。

exp[1:4,1:4]
##        GSM1366348 GSM1366349 GSM1366350 GSM1366351
## RFC2     8.932805   8.679543   8.625015   8.637085
## HSPA6    9.383421   8.605809   9.462774   9.898573
## PAX8     7.916751   8.500635   8.258467   8.553656
## GUCA1A   5.085221   2.414033   1.718570   4.311794

dim(exp)
## [1] 20161    22

kegg_list = split(KEGG_df$gene_symbol,KEGG_df$gs_exact_source)
lapply(kegg_list[1:3], head)

## $hsa00010
## [1] "ACSS1" "ACSS2" "ADH1A" "ADH1B" "ADH1C" "ADH4" 
## 
## $hsa00020
## [1] "ACLY" "ACO1" "ACO2" "CS"   "DLAT" "DLD" 
## 
## $hsa00030
## [1] "ALDOA" "ALDOB" "ALDOC" "DERA"  "FBP1"  "FBP2"

KEGG_ES <- gsva(expr=exp, 
               gset.idx.list=kegg_list, 
               parallel.sz=32) #自己电脑parallel.sz写5就好，线程数

## Setting parallel calculations through a MulticoreParam back-end
## with workers=32 and tasks=100.
## Estimating GSVA scores for 186 gene sets.
## Estimating ECDFs with Gaussian kernels
## Estimating ECDFs in parallel
## 
#   |=============================================| 100%

KEGG_ES[1:4,1:4]
##          GSM1366348 GSM1366349 GSM1366350  GSM1366351
## hsa00010 0.09619189 -0.2030580 -0.2575823  0.13014453
## hsa00020 0.41117392 -0.4726020 -0.4338775 -0.26959554
## hsa00030 0.03504000 -0.4179045 -0.4280300 -0.09500726
## hsa00040 0.16960156 -0.1391068 -0.2706679 -0.02557599

go_list = split(GO_df$gene_symbol,GO_df$gs_exact_source)
lapply(go_list[1:3], head)

## $`GO:0000002`
## [1] "AKT3"  "DNA2"  "FLCN"  "LIG3"  "LONP1" "MEF2A"
## 
## $`GO:0000003`
## [1] "AAAS"  "ABAT"  "ABCC2" "ABHD2" "ACE"   "ACOD1"
## 
## $`GO:0000012`
## [1] "AP002495.1" "APLF"       "APTX"       "ERCC6"     
## [5] "ERCC8"      "LIG4"

GO_ES <- gsva(expr=exp, 
               gset.idx.list=go_list, 
               parallel.sz=32) #自己电脑parallel.sz写5就好，线程数、

## Setting parallel calculations through a MulticoreParam back-end
## with workers=32 and tasks=100.
## Estimating GSVA scores for 10266 gene sets.
## Estimating ECDFs with Gaussian kernels
## Estimating ECDFs in parallel
## 
#   |=============================================| 100%

KEGG_ES[1:4,1:4]
##          GSM1366348 GSM1366349 GSM1366350  GSM1366351
## hsa00010 0.09619189 -0.2030580 -0.2575823  0.13014453
## hsa00020 0.41117392 -0.4726020 -0.4338775 -0.26959554
## hsa00030 0.03504000 -0.4179045 -0.4280300 -0.09500726
## hsa00040 0.16960156 -0.1391068 -0.2706679 -0.02557599

# 这里的参数kcdf，帮助文档里有说明：

# By default, kcdf="Gaussian" which is suitable when input expression values are continuous, such as microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input expression values are integer counts, such as those derived from RNA-seq experiments, then this argument should be set to kcdf="Poisson".

5.GSVA 差异分析
就是limma差异分析，敲简单

design = model.matrix(~Group)
fit = lmFit(GO_ES, design)
fit = eBayes(fit)
DEG = topTable(fit, coef = 2, number = Inf)
head(DEG)
##                 logFC     AveExpr          t
## GO:0070180 -1.0388245 -0.01574550 -10.988613
## GO:0000104 -0.9768231  0.01008804 -10.300079
## GO:0045273 -0.9768231  0.01008804 -10.300079
## GO:0010626 -0.8295852  0.01190610  -9.550725
## GO:0043426  0.8652926 -0.02684029   9.088616
## GO:0010499 -0.8009479  0.02692075  -9.065723
##                 P.Value    adj.P.Val        B
## GO:0070180 8.519107e-11 8.745715e-07 14.61388
## GO:0000104 3.066710e-10 1.049428e-06 13.42665
## GO:0045273 3.066710e-10 1.049428e-06 13.42665
## GO:0010626 1.318927e-09 3.385026e-06 12.06086
## GO:0043426 3.358846e-09 6.023685e-06 11.17890
## GO:0010499 3.520564e-09 6.023685e-06 11.13440

想做火山图、热图、PCA图，也都是可以借我的小函数一用：

draw_heatmap(GO_ES[head(rownames(DEG),200),],Group)
draw_volcano(DEG,pkg = 4,logFC_cutoff = 0.5)


draw_pca(GO_ES,Group)

