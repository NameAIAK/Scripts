#manifest
sample-id	forward-absolute-filepath	reverse-absolute-filepath
SRR33413006	/data/16stest/SRR33413006_1.fastq	/data/16stest/SRR33413006_2.fastq
SRR33413007	/data/16stest/SRR33413007_1.fastq	/data/16stest/SRR33413007_2.fastq
SRR33413008	/data/16stest/SRR33413008_1.fastq	/data/16stest/SRR33413008_2.fastq
SRR33413009	/data/16stest/SRR33413009_1.fastq	/data/16stest/SRR33413009_2.fastq


# meta
SampleID	group
SRR33413006	test1
SRR33413007	test1
SRR33413008	test2
SRR33413009	test2

# 下载数据
./prefetch SRR33413006
./prefetch SRR33413007

./prefetch SRR33413008
./prefetch SRR33413009

# 将数据拆分成fastq文件
fastq-dump --split-3 SRR33413006.sra


# 1
fastqc -f fastq -o
# 2
qiime tools import  --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest  --output-path paired-end-demux.qza  --input-format PairedEndFastqManifestPhred33V2
# 3
# qiime demux summarize --i-data paired-end-demux.qza  --o-visualization demux-summary.qzv
# 
qiime cutadapt trim-paired  --i-demultiplexed-sequences paired-end-demux.qza  --p-cores 16  --p-no-indels --p-front-f ACTCCTACGGGAGGCAGCAG --p-front-r GGACTACHVGGGTWTCTAAT --o-trimmed-sequences primer-trimmed-demux.qza

# qiime demux summarize --i-data primer-trimmed-demux.qza  --o-visualization primer-trimmed-demux.qzv

qiime vsearch merge-pairs  --i-demultiplexed-seqs primer-trimmed-demux.qza  --o-unmerged-sequences unmerged_demux-joined.qza --p-threads  16  --o-merged-sequences demux-joined.qza

# qiime demux summarize  --i-data demux-joined.qza  --o-visualization demux-joined-summary.qzv

qiime quality-filter q-score --p-min-quality 20 --i-demux demux-joined.qza --o-filtered-sequences demux-joined-filtered.qza --o-filter-stats demux-joined-filter-stats.qza

qiime metadata tabulate --m-input-file demux-joined-filter-stats.qza --o-visualization demux-joined-filter-stats.qzv

# qiime demux summarize --i-data demux-joined-filtered.qza --o-visualization demux-joined-filtered-summary.qzv

qiime deblur denoise-16S  --i-demultiplexed-seqs demux-joined-filtered.qza --p-trim-length 400 --p-sample-stats --p-jobs-to-start 2 --p-min-reads 1 --o-representative-sequences repset-seqs.qza --o-table feature-table.qza --o-stats deblur-stats.qza

# qiime feature-table summarize --i-table feature-table.qza --o-visualization feature-table.qzv

# qiime deblur visualize-stats --i-deblur-stats deblur-stats.qza --o-visualization deblur-stats.qzv                            

# 获取OTU表格：
# qiime feature-table summarize --i-table feature-table.qza --o-visualization  table.qzv --m-sample-metadata-file meta

qiime tools export --input-path feature-table.qza --output-path  dada2-table

biom convert -i feature-table.biom -o dada2-table/asv_dada2-table.txt  --table-type "OTU table" --to-tsv