import os
import argparse


parser = argparse.ArgumentParser(description="")
parser.add_argument('--i', type=str, help='map_reuslt')
args = parser.parse_args()
map_file=args.i

sample_name=map_file.split('.')[0]
map_reads = int(os.popen('cat '+map_file+' |cut -f1 |grep "C"|wc -l').read().strip())

r1='/data/metagenomic/rawdata/20240827/{0}/{1}.1.fastq.gz'.format(sample_name,sample_name)
r2='/data/metagenomic/rawdata/20240827/{0}/{1}.2.fastq.gz'.format(sample_name,sample_name)
raw_reads=int(int(os.popen('zcat ' + r1 + '/*gz |wc -l').read().strip())/4)

out=open('/data/metagenomic/workspace/20240828/map.stats','a')
ratio=map_reads/raw_reads
out.write(sample_name+'\t'+str(map_reads)+'\t'+str(raw_reads)+'\t'+str(ratio)+'\n')


