1. 质控

采用FASTQC查看测序数据质量

#fastqc -o FASTQC/ -t 8 Control_R1.fastq.gz Control_R2.fastq.gz Treated_R1.fastq.gz Treated _R2.fastq.gz

#multiqc ./



2. 过滤

采用Cutadapt对测序文件进行过滤，目的包括：去除测序引物及接头、去除reads两端低质量碱基、去除N碱基过多的reads、去除截短后单端reads长度小于75bp的reads。

#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 30 -m 75 --trim-n --report=minimal -o Control_out_R1.fastq.gz -p Control_out_R2.fastq.gz Control_R1.fastq.gz Control_R2.fastq.gz



3.比对及排序

使用Bowtie2将clean data与参考基因组进行比对并采用samtools对BAM文件进行排序。其中-X代表最大插入片段，宽泛的插入范围为10-1000bp，一般一个核小体147 180，大片段很稀有。

#bowtie2 -p 10 -X 1000 -x bowtie_index -1 Control_out_R1.fastq -2 Control_out_R2.fastq |samtools sort -O bam -@ 5 -o Control.bam



4. 采用sambamba及samtools去除PCR重复以及线粒体基因，去除低质量序列

#sambamba markdup -r Control.bam Control.sambamba.rmdup.bam 去除PCR重复

#samtools view -h -f 2 -q 30 Control.sambamba.rmdup.bam |grep -v chrM |samtools sort -O bam -@ 5 -o -> Control.last.bam 去除线粒体重复，去除低质量序列



5. samtools对测序深度、覆盖度、比对率、重复率等进行统计

#samtools index Control.sambamba.rmdup.bam

#samtools flagstat Control.sambamba.rmdup.bam > Control.rmdup.stat

#cat Control.rmdup.stat 查看stat内容



6. 采用bedtools将bam文件转换为bed文件

#bedtools bamtobed -i Control.last.bam >Control.last.bed

#less Control.last.bed 查看bed文件

#bedtools intersect -a ../peaks/Control.bed -b Control_summits.bed 计算插入片段长度



7.采用macs2进行call peaks

# macs2 callpeak -t Control.last.bed -g hs --nomodel --shift -100 --extsize 200 -n Control --outdir ../peaks

# wc Control_peaks.narrowPeak 查看样本peaks

# wc -l *bed 查看所有样本的peaks数



8. 采用deeptools及IGV进行可视化

对bam文件进行归一化之后可进行IGV可视化，同时将bam文件转化为bw文件

#ls *last.bam |while read id; do

> nohup bamCoverage --normalizeUsing CPM -b $id -o ${id%%.*}.last.bw &

> done |




#mm data    SRR14305716 SRR14305717	 
# 下载数据
../16stest/sratoolkit.3.0.2-ubuntu64/bin/prefetch SRR14305716
../16stest/sratoolkit.3.0.2-ubuntu64/bin/prefetch SRR14305717

# analysis  https://zhuanlan.zhihu.com/p/415718382
# QC
less atac_ara_test.txt |while  read id;
do 
fastqc ${id}_1.fastq -o /path/fastqc
fastqc ${id}_2.fastq -o /path/fastqc
done

less /home/xuyu/atac_ara/atac_ara_test.txt |while  read id;
do 
fastp -i ${id}_1.fastq -o ./fastp/${id}_1.fastq -I ${id}_2.fastq -O ./fastp/${id}_2.fastq -f 16 -t 2 -L
done

## conda装包
conda install -c bioconda bwa
## 建基因组index
bwa index -a bwtsw /home/xuyu/atac_ara/atac_data/genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa



less /home/xuyu/atac_ara/atac_ara_test.txt |while  read id;
do
bwa mem -v 3 -t 4  /home/xuyu/atac_ara/atac_data/genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa   ./${id}_1.fastq ./${id}_2.fastq -o ./bwa_sam/${id}_bwa.sam
done

# sam2bam
less /home/xuyu/atac_ara/atac_ara_test.txt |while  read id;
do
samtools view -@ 4 ./bwa_sam/${id}_bwa.sam -bF 12 -q 10  -O bam -o ./bwa_bam_samtools10/${id}_samtools10.bam
done

# 去重复
less /home/xuyu/atac_ara/atac_ara_test.txt |while  read id;
do
sambamba markdup -r -t 4 ${id}_samtools10.bam ${id}_samtools10_rdup.bam
done













