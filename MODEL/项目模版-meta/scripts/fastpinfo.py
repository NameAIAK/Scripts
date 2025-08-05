import os
import json

jsonfolder = './fastp_json'

fw= open('fastp_summary.txt', 'a')
fw.write('samplename\tafter_filter_readsT\tq20_rate\tq30_rate\tgc_content\tafter_filter_reads1\tafter_filter_reads2\tread1_mean_length\tread2_mean_length\n')

for root,dirs,files in os.walk(jsonfolder):
    for file in files:
        if file.endswith('.json'):
            jsonfile = os.path.join(root, file)
            print(jsonfile)
            samplename = file.split('.json')[0]
            with open(jsonfile, 'r') as f:
                data = json.load(f)
                after_filter_readsT=str(data['summary']['after_filtering']['total_reads'])
                q20_rate=str(data['summary']['after_filtering']['q20_rate'])
                q30_rate=str(data['summary']['after_filtering']['q30_rate'])
                read1_mean_length=str(data['summary']['after_filtering']['read1_mean_length'])
                read2_mean_length=str(data['summary']['after_filtering']['read2_mean_length'])
                gc_content=str(data['summary']['after_filtering']['gc_content'])
                after_filter_reads1=str(data['read1_after_filtering']['total_reads'])
                after_filter_reads2=str(data['read2_after_filtering']['total_reads'])
        info=samplename+'\t'+after_filter_readsT+'\t'+q20_rate+'\t'+q30_rate+'\t'+gc_content+'\t'+after_filter_reads1+'\t'+after_filter_reads2+'\t'+read1_mean_length+'\t'+read2_mean_length+'\n'
        fw.write(info)