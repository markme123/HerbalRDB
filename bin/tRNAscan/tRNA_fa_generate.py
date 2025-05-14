# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 13:47:34 2019

@author: Administrator
"""

import os
import argparse
import re

parser = argparse.ArgumentParser(description='generate tRNA.fa')
parser.add_argument('-gff', dest='gff', type=str,
                    help='tRNA gff')
parser.add_argument('-fa', dest='fa', type=str,
                    help='genome fa', required=True)
args = parser.parse_args()



fa_dict = {}
trna_dict = {}
chr_name = []
chr_seq = []
tmp_line = ""
with open(args.fa) as f2:
    for line in f2:
        if line.startswith('>'):
            chr_name.append(line.strip()[1:].split(' ')[0])
            if tmp_line != "":
                chr_seq.append(tmp_line)
            tmp_line =""
        else:
            tmp_line += line.strip()
    chr_seq.append(tmp_line)
for i in range(len(chr_name)) :
    fa_dict[chr_name[i]] = chr_seq[i]           
            
            
with open(args.gff) as f1:
    f1.readline()
    for line in f1:
        getid = re.search('ID=(.*?);',line)
        if getid:
            name = getid.group(1)
            tmp_list = line.strip().split('\t')
            trna_dict[tmp_list[0]+'.'+name] = (fa_dict[tmp_list[0]])[int(tmp_list[3])-1:int(tmp_list[4])]

with open(os.path.basename(args.fa)+'tRNAdb','w') as f:
    for key in trna_dict.keys():
        f.write('>'+key+'\n')
        f.write(trna_dict[key]+'\n')

        
