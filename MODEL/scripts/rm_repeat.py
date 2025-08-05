#!/usr/bin/env python3
# -*- encoding: utf-8 -*-



import argparse
from Bio import SeqIO
import time
# python rm_repeat.py -i rawCRISPRS.fna -r repeat.txt -o no_repeat.fasta

parser = argparse.ArgumentParser(prog="",
                                     description="")

parser.add_argument("-i", "--rawCRISPRS", dest="",
                    required=True, type=str, help="")
parser.add_argument("-r", "--repeatfile", dest="",
                    required=True, type=str, help="")
parser.add_argument("-o", "--CRISPRS", dest="", 
                    type=str, help="")
args = parser.parse_args()


rawCRISPRS=args.rawCRISPRS
repeatfile=args.repeatfile
outCRISPRS=args.CRISPRS

repeat_list=[]
with open(repeatfile) as f:
    for line in f:
        line = line.strip()
        repeat_list.append(line)
# print(repeat_list)

records = SeqIO.parse(rawCRISPRS, "fasta")
corrected = open(outCRISPRS, 'a')
for record in records:
    for repeat in repeat_list:
        if repeat in record.seq:
            # print(record.seq)
            record.seq=(record.seq).replace(repeat,'')
            # print(record.seq)
            
            # time.sleep(3)
            break
        else:
            continue
    SeqIO.write(record, corrected, "fasta")