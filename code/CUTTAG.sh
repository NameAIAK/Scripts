# conda activate base
# # 新建分析结果文件夹
# mkdir 1_rawdata 2_clean 3_sam 3_bam 4_rmdup 5_bw 6_callpeak
# #  复制数据到1_rawdata；数据命名：N87_1_R1.fq.gz    N87_1_R2.fq.gz  sample_R1.fq.gz sample_R2.fq.gz
# fileendr1='_clean_R1.fq.gz'
# fileendr2='_clean_R2.fq.gz'
# ls *$fileendr1 |while read id ; do cp $id $projectname/1_rawdata/${id%$fileendr1}_R1.fq.gz; done
# ls *$fileendr2 |while read id ; do cp $id $projectname/1_rawdata/${id%$fileendr2}_R2.fq.gz; done


#!/bin/bash
echo "Finish Copy & Prepare..."
# 脚本使用方法:bash CUTTAG.sh /data/CUTTag/n87 human N87）
# 进入到项目分析文件夹
projectname=$1
cd $projectname
#human/mouse
datatype=$2
pre=$3

# 1. clean data : 
conda activate RNA-seq
ls 1_rawdata/*_R1.fq.gz |while read id ; do trim_galore -j 50 -q 30 --phred33 --length 120 -e 0.1  --paired -o 2_clean/ $id ${id%_R1.*gz}_R2.fq.gz;done

# 2.map: 
#小鼠数据库/data/RNA/reference/bowtie2_mouse/mm10；人源数据库/data/RNA/reference/bowtie2_human/GRCh38/GRCh38
cd 2_clean

if [[ "$datatype" == "mouse" ]]; then
    echo "Running bowtie2 for MOUSE (mm10 reference)..."
    ls *R1*.fq.gz  | while read id ; do bowtie2  -p 50  -x /data/RNA/reference/bowtie2_mouse/mm10  --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -1 ${id%_R1*}_R1_val_1.fq.gz -2 ${id%_R1*}_R2_val_2.fq.gz  -S ../3_sam/${id%_R1*}.sam;done
elif [ "$datatype" == "human" ]; then
    echo "Running bowtie2 for HUMAN (GRCh38 reference)..."
    #人人人人人人人人人人人人人人人人人人人人
    ls *R1*.fq.gz  | while read id ; do bowtie2  -p 50  -x /data/RNA/reference/bowtie2_human/GRCh38/GRCh38  --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -1 ${id%_R1*}_R1_val_1.fq.gz -2 ${id%_R1*}_R2_val_2.fq.gz  -S ../3_sam/${id%_R1*}.sam;done
fi

# 3.samtobam:
cd ../3_sam
ls *.sam |while read id ;do samtools view -@ 50 -bF 12 -q 30 -S $id > ../3_bam/${id%.sam}.map.q30.F12.bam;done

# 3-2.remove duplication: 
cd ../3_bam
ls *.bam |while read id ;do /home/star/anaconda3/envs/sambamba/bin/sambamba markdup -r $id ../4_rmdup/${id%.bam}.rmdup.bam ;done

# #计算去除重复后的序列数
# samtools view N87_1.map.q30.F12.rmdup.bam |less -S|wc -l 

# 4 sort bam and make index
cd ../4_rmdup
ls *.rmdup.bam |while read id ;do samtools sort -@ 50 $id -o ${id%.rmdup.bam}.rmdup.sorted.bam;done
ls *.rmdup.sorted.bam |while read id ;do samtools index $id;done

# 5.bam to bw: 
ls *.rmdup.sorted.bam |while read id ;do  bamCoverage -b $id --normalizeUsing RPKM -p 50 -o ../5_bw/${id%.map.q30.F12.rmdup.sorted.bam}.bw ; done

# 6.call peak(adjust for transcrip factor/ Histone modification type)：
conda activate macs2
macs2  callpeak -t *.map.q30.F12.rmdup.sorted.bam --bdg -p 1e-5 -g hs -n '../6_callpeak/'$pre'.peak'

# 7.visualization
conda activate ChIP-seq
cd ../6_callpeak
computeMatrix  reference-point -p 20 -R $pre'.peak_summits.bed' -a 3000 -b 3000 -S  ../5_bw/*.bw --skipZeros  -out ./$pre'.computeMatrix.gz'
#通过此图可以查看样本测序情况
plotHeatmap -m $pre'.computeMatrix.gz' -o $pre'.computeMatrix.gz.pdf' --colorMap RdBu --zMin -3 --zMax 3

# 8.annotation
conda activate macs2
annotatePeaks.pl $pre'.peak_summits.bed' hg38 > $pre'_peak.anno.txt'
tail +2 '../6_callpeak/'$pre'_peak.anno.txt' |cut -f 2,3,4 > '../6_callpeak/'$pre'_peak.txt'

# # 9.bam2tdf
conda activate RNA-seq
ls ../4_rmdup/*.map.q30.F12.rmdup.sorted.bam | while read id ;do igvtools count $id ${id%.map.q30.F12.rmdup.sorted.bam}'.tdf' hg38.chrom.sizes 1>>bam2tdf.log 2>&1;done
# igvtools count ../4_rmdup/$sample1'.map.q30.F12.rmdup.sorted.bam '$pre'_1.tdf' hg38.chrom.sizes 1>>bam2tdf.log 2>&1
# igvtools count ../4_rmdup/$sample2'.map.q30.F12.rmdup.sorted.bam '$pre'_2.tdf' hg38.chrom.sizes 1>>bam2tdf.log 2>&1

##############################个性化分析##############################
# 合并bam文件
# samtools merge -o AGS_merge.bam ../4_rmdup/AGS_1_.map.q30.F12.rmdup.sorted.bam ../4_rmdup/AGS_2_.map.q30.F12.rmdup.sorted.bam
# 生成索引
# samtools index AGS_merge.bam
# 生成bw文件
# bamCoverage -b ../AGS_merge.bam --normalizeUsing RPKM -p 50 -o AGS_M.bw

# conda activate macs2
# macs2  callpeak -t *bam --bdg -p 1e-5 -g hs -n AGS_N87_M.peak
# conda activate ChIP-seq
# computeMatrix  reference-point -p 20 -R AGS_N87_M.peak_summits.bed -a 3000 -b 3000 -S *.bw --skipZeros  -out AGS_N87_M.computeMatrix.gz
# #通过此图可以查看样本测序情况
# plotHeatmap -m AGS_N87_M.computeMatrix.gz -o AGS_N87_M.computeMatrix.gz.pdf --colorMap RdBu --zMin -3 --zMax 3

# # inner
# conda activate RNA-seq
# bedtools intersect -a AGS.peak_summits.bed -b N87.peak_summits.bed > AGS_N87_overlaps.bed
# sort -k1,1 -k2,2n AGS_N87_overlaps.bed > AGS_N87_overlaps.sorted.bed
# bedtools merge -i AGS_N87_overlaps.sorted.bed > AGS_N87_overlaps.sorted.merge.bed

# conda activate macs2
# annotatePeaks.pl AGS_N87_M.peak_summits.bed hg38 > AGS_N87_M_peak.anno.txt

# conda activate RNA-seq
# ls *_M.bam | while read id ;do igvtools count $id ${id%.bam}'.tdf' hg38.chrom.sizes 1>>bam2tdf.log 2>&1;done