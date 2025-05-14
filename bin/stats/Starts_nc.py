# -*- coding: utf-8 -*-
"""
Created on Thu May 16 16:28:32 2019

@author: Administrator
"""

import os
import re
import argparse
import gzip

parser = argparse.ArgumentParser(description='Starts gff')
parser.add_argument('-m',dest='miRNA',type=str,help='miRNA gff',default='./miRNA.gff')
parser.add_argument('-t',dest='tRNA',type=str,help='tRNA gff',default='./tRNA.gff')
parser.add_argument('-r',dest='rRNA',type=str,help='rRNA gff',default='./rRNA.gff')
parser.add_argument('-s',dest='snRNA',type=str,help='snRNA gff',default='./snRNA.gff')
parser.add_argument('-g',dest='ref',type=str,help='genome.fasta',default='./genome.fasta')
parser.add_argument('-o',dest='output',type=str,help='output dir',default='./')
parser.add_argument('-sp',dest='species',type=str,help='species name',default='anno')
args = parser.parse_args()

def genome_size(ref):
    genome_size = 0
    if not os.path.exists(ref):
        print('Please input right genome file !!!!')
        exit(1)
    if ref.endswith('.gz'):
        with gzip.open(ref) as f:
            for line in f:
                if line.startswith('>'):
                    pass
                else:
                    genome_size += len(line.strip())
    else:
         with open(ref) as f:
            for line in f:
                if line.startswith('>'):
                    pass
                else:
                    genome_size += len(line.strip())
    if genome_size == 0:
        print('Please input right genome file !!!!')
        exit(1)
    else:
        return genome_size

total_size = genome_size(args.ref)
miRNA_dict = {}
miRNA_dict['miRNA'] = {}
miRNA_dict['miRNA']['Copy'] = 0
miRNA_dict['miRNA']['Average'] = 0
miRNA_dict['miRNA']['Total'] = 0
miRNA_dict['miRNA']['genome'] = 0
with open(args.miRNA,'r') as f1:
    for line in f1:
        if line.strip() == '' or line.startswith('#'):
            pass
        else:
            miRNA_dict['miRNA']['Copy'] += 1
            len_m = int(line.split('\t')[4]) - int(line.split('\t')[3]) + 1
            miRNA_dict['miRNA']['Total'] += len_m
if miRNA_dict['miRNA']['Copy'] == 0:
    miRNA_dict['miRNA']['Average'] = 0
else:
    miRNA_dict['miRNA']['Average'] = miRNA_dict['miRNA']['Total']/miRNA_dict['miRNA']['Copy']
#    miRNA_dict['miRNA']['Average'] = format(float(miRNA_dict['miRNA']['Total'])/float(miRNA_dict['miRNA']['Copy']),'.4f')
miRNA_dict['miRNA']['genome'] = miRNA_dict['miRNA']['Total']*100/total_size
tRNA_dict = {}
tRNA_dict['tRNA'] = {}
tRNA_dict['tRNA']['Copy'] = 0
tRNA_dict['tRNA']['Average'] = 0
tRNA_dict['tRNA']['Total'] = 0
tRNA_dict['tRNA']['genome'] = 0
with open(args.tRNA,'r') as f1:
    for line in f1:
        if line.strip() == '' or line.startswith('#'):
            pass
        else:
            if 'Type=Pseudo' in line.split('\t')[-1]:
                pass
            else:
                tRNA_dict['tRNA']['Copy'] += 1
                len_m = int(line.split('\t')[4]) - int(line.split('\t')[3]) + 1
                tRNA_dict['tRNA']['Total'] += len_m
if tRNA_dict['tRNA']['Copy'] == 0:
    tRNA_dict['tRNA']['Average'] = 0
else:
    tRNA_dict['tRNA']['Average'] = tRNA_dict['tRNA']['Total']/tRNA_dict['tRNA']['Copy']
#    tRNA_dict['tRNA']['Average'] = format(float(tRNA_dict['tRNA']['Total'])/float(tRNA_dict['tRNA']['Copy']),'.4f')
tRNA_dict['tRNA']['genome'] = tRNA_dict['tRNA']['Total']*100/total_size

rRNA_dict = {}
rRNA_dict['rRNA'] = {}
rRNA_dict['rRNA']['Copy'] = 0
rRNA_dict['rRNA']['Average'] = 0
rRNA_dict['rRNA']['Total'] = 0
rRNA_dict['rRNA']['genome'] = 0
#rRNA_list = ['18S','28S','5.8S','5S','8S']
rRNA_list = ['18S','28S','5.8S','5S']
for i in rRNA_list:
    rRNA_dict[i] = {}
    rRNA_dict[i]['Copy'] = 0
    rRNA_dict[i]['Average'] = 0
    rRNA_dict[i]['Total'] = 0
    rRNA_dict[i]['genome'] = 0
with open(args.rRNA,'r') as f1:
    for line in f1:
        if line.strip() == '' or line.startswith('#'):
            pass
        else:
            rRNA_dict['rRNA']['Copy'] += 1
            len_m = int(line.split('\t')[4]) - int(line.split('\t')[3]) + 1
            rRNA_dict['rRNA']['Total'] += len_m
            for j in rRNA_list:
                if j in line.split('\t')[-1]:
                    rRNA_dict[j]['Copy'] += 1
                    rRNA_dict[j]['Total'] += len_m
                    break
                else:
                    pass
if rRNA_dict['rRNA']['Copy'] == 0:
    rRNA_dict['rRNA']['Average'] = 0
else:
    rRNA_dict['rRNA']['Average'] = rRNA_dict['rRNA']['Total']/rRNA_dict['rRNA']['Copy']
#    rRNA_dict['rRNA']['Average'] = format(float(rRNA_dict['rRNA']['Total'])/float(rRNA_dict['rRNA']['Copy']),'.4f')
rRNA_dict['rRNA']['genome'] = rRNA_dict['rRNA']['Total']*100/total_size
for i in rRNA_list:
    if rRNA_dict[i]['Copy'] == 0:
        rRNA_dict[i]['Average'] = 0
    else:
        rRNA_dict[i]['Average'] = rRNA_dict[i]['Total']/rRNA_dict[i]['Copy']
#        rRNA_dict[i]['Average'] = format(float(rRNA_dict[i]['Total'])/float(rRNA_dict[i]['Copy']),'.4f')
    rRNA_dict[i]['genome'] = rRNA_dict[i]['Total']*100/total_size

snRNA_dict = {}
snRNA_dict['snRNA'] = {}
snRNA_dict['snRNA']['Copy'] = 0
snRNA_dict['snRNA']['Average'] = 0
snRNA_dict['snRNA']['Total'] = 0
snRNA_dict['snRNA']['genome'] = 0
snRNA_list = ['CD-box','HACA-box','splicing','scaRNA']
for i in snRNA_list:
    snRNA_dict[i] = {}
    snRNA_dict[i]['Copy'] = 0
    snRNA_dict[i]['Average'] = 0
    snRNA_dict[i]['Total'] = 0
    snRNA_dict[i]['genome'] = 0
with open(args.snRNA,'r') as f1:
    for line in f1:
        if line.strip() == '' or line.startswith('#'):
            pass
        else:
            snRNA_dict['snRNA']['Copy'] += 1
            len_m = int(line.split('\t')[4]) - int(line.split('\t')[3]) + 1
            snRNA_dict['snRNA']['Total'] += len_m
            for j in snRNA_list:
                if j in line.split('\t')[-1]:
                    snRNA_dict[j]['Copy'] += 1
                    snRNA_dict[j]['Total'] += len_m
                    break
                else:
                    pass
if snRNA_dict['snRNA']['Copy'] == 0:
    snRNA_dict['snRNA']['Average'] = 0
else:
    snRNA_dict['snRNA']['Average'] = snRNA_dict['snRNA']['Total']/snRNA_dict['snRNA']['Copy']
#    snRNA_dict['snRNA']['Average'] = format(float(snRNA_dict['snRNA']['Total'])/float(snRNA_dict['snRNA']['Copy']),'.4f')
snRNA_dict['snRNA']['genome'] = snRNA_dict['snRNA']['Total']*100/total_size
for j in snRNA_list:
    if snRNA_dict[j]['Copy'] == 0:
        snRNA_dict[j]['Average'] = 0
    else:
        snRNA_dict[j]['Average'] = snRNA_dict[j]['Total']/snRNA_dict[j]['Copy']
#        snRNA_dict[j]['Average'] = format(float(snRNA_dict[j]['Total'])/float(snRNA_dict[j]['Copy']),'.4f')
    snRNA_dict[j]['genome'] = snRNA_dict[j]['Total']*100/total_size

nc_starts = os.path.join(args.output,args.species + '.ncRNA.statistics.xls')
with open(nc_starts,'w') as fw:
#    fw.write('Total Gene Lenght = %s\n\n'%(total_size))
    fw.write('Type\t\tCopy\tAverage length(bp)\tTotal length(bp)\t% of genome\n')
    fw.write('miRNA\t\t%s\t%.0f\t%s\t%.4f\n'%(miRNA_dict['miRNA']['Copy'],miRNA_dict['miRNA']['Average'],miRNA_dict['miRNA']['Total'],miRNA_dict['miRNA']['genome']))
    fw.write('tRNA\t\t%s\t%.0f\t%s\t%.4f\n'%(tRNA_dict['tRNA']['Copy'],tRNA_dict['tRNA']['Average'],tRNA_dict['tRNA']['Total'],tRNA_dict['tRNA']['genome']))
    fw.write('rRNA\trRNA\t%s\t%.0f\t%s\t%.4f\n'%(rRNA_dict['rRNA']['Copy'],rRNA_dict['rRNA']['Average'],rRNA_dict['rRNA']['Total'],rRNA_dict['rRNA']['genome']))
    for i in rRNA_list:
        fw.write('\t%s\t%s\t%.0f\t%s\t%.4f\n'%(i,rRNA_dict[i]['Copy'],rRNA_dict[i]['Average'],rRNA_dict[i]['Total'],rRNA_dict[i]['genome']))
    fw.write('snRNA\tsnRNA\t%s\t%.0f\t%s\t%.4f\n'%(snRNA_dict['snRNA']['Copy'],snRNA_dict['snRNA']['Average'],snRNA_dict['snRNA']['Total'],snRNA_dict['snRNA']['genome']))
    for j in snRNA_list:
        fw.write('\t%s\t%s\t%.0f\t%s\t%.4f\n'%(j,snRNA_dict[j]['Copy'],snRNA_dict[j]['Average'],snRNA_dict[j]['Total'],snRNA_dict[j]['genome']))
tmp_nc_starts = os.path.join(args.output,args.species + '.ncRNA.statistics_NoType.xls')
with open(tmp_nc_starts,'w') as fw:
    fw.write('Type\tCopy\tAverage length(bp)\tTotal length(bp)\t% of genome\n')
    fw.write('miRNA\t\t%s\t%.0f\t%s\t%.4f\n'%(miRNA_dict['miRNA']['Copy'],miRNA_dict['miRNA']['Average'],miRNA_dict['miRNA']['Total'],miRNA_dict['miRNA']['genome']))
    fw.write('tRNA\t%s\t%.0f\t%s\t%.4f\n'%(tRNA_dict['tRNA']['Copy'],tRNA_dict['tRNA']['Average'],tRNA_dict['tRNA']['Total'],tRNA_dict['tRNA']['genome']))
    fw.write('rRNA\t%s\t%.0f\t%s\t%.4f\n'%(rRNA_dict['rRNA']['Copy'],rRNA_dict['rRNA']['Average'],rRNA_dict['rRNA']['Total'],rRNA_dict['rRNA']['genome']))
    fw.write('snRNA\t%s\t%.0f\t%s\t%.4f\n'%(snRNA_dict['snRNA']['Copy'],snRNA_dict['snRNA']['Average'],snRNA_dict['snRNA']['Total'],snRNA_dict['snRNA']['genome']))
