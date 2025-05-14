# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:15:13 2019

@author: Administrator
"""

import os
import re
import argparse

parser = argparse.ArgumentParser(description='Merge blast and cmscan result')
parser.add_argument('-icm',dest='incmscan',type=str,help='cmscan gff file',default='./rRNA_cmscan.gff')
parser.add_argument('-ibl',dest='inblast',type=str,help='blast gff file',default='./rRNA_blast.gff')
parser.add_argument('-o',dest='output',type=str,help='output dir',default='./')
parser.add_argument('-s',dest='species',type=str,help='species name',default='anno')
args = parser.parse_args()

def sort_gff(gff_list):
    chr_list = []
    chr_dict = {}
    chr_start_dict = {}
    num = {}
    num['snRNA'] = 0
    num['rRNA'] = 0
    num['miRNA'] = 0
    new_gff_list = []
    for i in gff_list:
        if i.split('\t')[0] not in chr_list:
            chr_list.append(i.split('\t')[0])
        else:
            pass
        if i.split('\t')[0] not in chr_dict:
            chr_start_dict[i.split('\t')[0]] = []
            chr_start_dict[i.split('\t')[0]].append(float(i.split('\t')[3] + '.' + i.split('\t')[4]))
            chr_dict[i.split('\t')[0]] = {}
            if float(i.split('\t')[3] + '.' + i.split('\t')[4]) not in chr_dict[i.split('\t')[0]]:
                chr_dict[i.split('\t')[0]][float(i.split('\t')[3] + '.' + i.split('\t')[4])] = []
                chr_dict[i.split('\t')[0]][float(i.split('\t')[3] + '.' + i.split('\t')[4])].append(i)
            else:
                chr_dict[i.split('\t')[0]][float(i.split('\t')[3] + '.' + i.split('\t')[4])].append(i)
        else:
            chr_start_dict[i.split('\t')[0]].append(float(i.split('\t')[3] + '.' + i.split('\t')[4]))
            if float(i.split('\t')[3] + '.' + i.split('\t')[4]) not in chr_dict[i.split('\t')[0]]:
                chr_dict[i.split('\t')[0]][float(i.split('\t')[3] + '.' + i.split('\t')[4])] = []
                chr_dict[i.split('\t')[0]][float(i.split('\t')[3] + '.' + i.split('\t')[4])].append(i)
            else:
                chr_dict[i.split('\t')[0]][float(i.split('\t')[3] + '.' + i.split('\t')[4])].append(i)            
    for j in chr_list:
        set_list = list(set(chr_start_dict[j]))
        set_list.sort()
        for k in set_list:
            for l in chr_dict[j][k]:
                num[l.split('\t')[2]] += 1
                new_gff_list.append('\t'.join(l.split('\t')[:-1]) + '\t' + 'ID=%s_%s;%s'%(l.split('\t')[2],str(num[l.split('\t')[2]]),';'.join(l.split('\t')[-1].split(';')[1:])))
    return new_gff_list

gff_file = []
scan_5s = []
scan_5_8s = []
flag = 0
tag = 0
if not os.path.exists(args.incmscan):
    print('Please input right cmscan_gff !')
else:
    with open(args.incmscan,'r') as fs:
        for line in fs:
            if line.startswith('#') or line.strip() == '':
                pass
            else:
                if '5S' in line.strip().split('\t')[-1].split('Annotation=')[1]:
                    scan_5s.append(line.strip().split('\t')[3] + '_' + line.strip().split('\t')[4])
                    gff_file.append(line.strip())
                    flag += 1
                elif '5.8S' in line.strip().split('\t')[-1].split('Annotation=')[1]:
                    scan_5_8s.append(line.strip().split('\t')[3] + '_' + line.strip().split('\t')[4])
                    gff_file.append(line.strip())
                    flag += 1
                else:
                    pass
if not os.path.exists(args.inblast):
    print('Please input right blast_gff !')
else:
    with open(args.inblast,'r') as fi:
        for line in fi:
            if line.startswith('#') or line.strip() == '':
                pass
            else:
                if '5S' in line.strip().split('\t')[-1].split('Annotation=')[1]:
                    for i in scan_5s:
                        if (int(line.strip().split('\t')[3]) >= int(i.split('_')[0]) and int(line.strip().split('\t')[3]) <= int(i.split('_')[1])) or (int(line.strip().split('\t')[4]) >= int(i.split('_')[0]) and int(line.strip().split('\t')[4]) <= int(i.split('_')[1])):
                            tag = 1
                        else:
                            pass
                    if tag == 1:
                        pass
                    else:
                        if 'ID=' in line.strip().split('\t')[-1]:
                            gff_file.append('\t'.join(line.strip().split('\t')[:-1]) + '\t' + 'ID=rRNA_' + str(flag) + ';' + ';'.join(line.strip().split('\t')[-1].split('ID=')[1].split(';')[1:]))
                        else:
                            gff_file.append('\t'.join(line.strip().split('\t')[:-1]) + '\t' + 'ID=rRNA_' + str(flag) + ';' + line.strip().split('\t')[-1])
                elif '5.8S' in line.strip().split('\t')[-1].split('Annotation=')[1]:
                    for i in scan_5_8s:    
                        if (int(line.strip().split('\t')[3]) >= int(i.split('_')[0]) and int(line.strip().split('\t')[3]) <= int(i.split('_')[1])) or (int(line.strip().split('\t')[4]) >= int(i.split('_')[0]) and int(line.strip().split('\t')[4]) <= int(i.split('_')[1])):
                            tag = 1
                        else:
                            pass
                    if tag == 1:
                        pass
                    else:
                        if 'ID=' in line.strip().split('\t')[-1]:
                            gff_file.append('\t'.join(line.strip().split('\t')[:-1]) + '\t' + 'ID=rRNA_' + str(flag) + ';' + ';'.join(line.strip().split('\t')[-1].split('ID=')[1].split(';')[1:]))
                        else:
                            gff_file.append('\t'.join(line.strip().split('\t')[:-1]) + '\t' + 'ID=rRNA_' + str(flag) + ';' + line.strip().split('\t')[-1])
                else:
                    flag += 1
                    if 'ID=' in line.strip().split('\t')[-1]:
                        gff_file.append('\t'.join(line.strip().split('\t')[:-1]) + '\t' + 'ID=rRNA_' + str(flag) + ';' + ';'.join(line.strip().split('\t')[-1].split('ID=')[1].split(';')[1:]))
                    else:
                        gff_file.append('\t'.join(line.strip().split('\t')[:-1]) + '\t' + 'ID=rRNA_' + str(flag) + ';' + line.strip().split('\t')[-1])
if not os.path.exists(args.output):
    os.makedirs(args.output)
    args.output = os.patn.abspath(args.output)
wfile = os.path.join(args.output,args.species + '.rRNA.gff')
with open(wfile,'w') as fw:
    fw.write('##gff-version 3\n')
    fw.write('\n'.join(sort_gff(gff_file)) + '\n')