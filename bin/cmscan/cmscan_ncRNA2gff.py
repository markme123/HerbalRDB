# -*- coding: utf-8 -*-
"""
Created on Thu May 16 16:28:32 2019

@author: Administrator
"""

import os
import re
import argparse

parser = argparse.ArgumentParser(description='Change cmscan result to gff')
#parser.add_argument('-nctype',dest='nctype',type=str,help='ncRNA type(rRNA,miRNA,rRNA,snRNA)',default='all')
parser.add_argument('-i',dest='input',type=str,help='cmscan result file',default='./ga.tblout')
parser.add_argument('-o',dest='output',type=str,help='output dir',default='./')
parser.add_argument('-t',dest='Rfam_anno',type=str,help='Rfam_anno.txt',default='/home/soft_env/Rfam/newRfam')
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
type_nc = {}
Annotation = {}
if not os.path.exists(args.Rfam_anno):
    print('Please input right Rfam_anno.txt !')
else:
    with open(args.Rfam_anno,'r') as fm:
        for line in fm:
            if line.strip() == '' or line.startswith('#'):
                pass
            else:
                l = line.strip().split('\t')
                type_nc[l[0]] = l[2]
                Annotation[l[0]] = l[3]
                
gff_dict = {}
gff_dict['snRNA'] = []
gff_dict['rRNA'] = []
gff_dict['miRNA'] = []
num = {}
num['snRNA'] = 0
num['rRNA'] = 0
num['miRNA'] = 0
if not os.path.exists(args.Rfam_anno):
    print('Please input right Rfam_anno.txt !')
else:
    with open(args.input,'r') as fr:
        for line in fr:
            if line.startswith('#'):
                pass
            elif line.strip() == '':
                pass
            else:
                gff_line = []
                l = line.strip()
                desc = re.sub(r'\s+','\t',l)
                perm = desc.split('\t')
                if perm[2] in type_nc:
                    type_line = type_nc[perm[2]].split('; ')[1].replace(';','')
                    num[type_line] += 1
                    if perm[11] == '-':
                        gff_line.append(perm[3])
                        gff_line.append('cmscan')
                        gff_line.append(type_line)
                        gff_line.append(perm[10])
                        gff_line.append(perm[9])
                        gff_line.append(perm[17])
                        gff_line.append(perm[11])
                        gff_line.append('.')
                        gff_line.append('ID=%s_%s;Target=%s %s %s;type="%s";Annotation=%s'%(type_line,str(num[type_line]),perm[2],perm[7],perm[8],type_nc[perm[2]],Annotation[perm[2]]))
                    else:
                        gff_line.append(perm[3])
                        gff_line.append('cmscan')
                        gff_line.append(type_line)
                        gff_line.append(perm[9])
                        gff_line.append(perm[10])
                        gff_line.append(perm[17])
                        gff_line.append(perm[11])
                        gff_line.append('.')
                        gff_line.append('ID=%s_%s;Target=%s %s %s;type="%s";Annotation=%s'%(type_line,str(num[type_line]),perm[2],perm[7],perm[8],type_nc[perm[2]],Annotation[perm[2]]))
                    gff_dict[type_line].append('\t'.join(gff_line))
                else:
                    pass
if not os.path.exists(args.output):
    os.makedirs(args.output)
    args.output = os.patn.abspath(args.output)
for nc in gff_dict:
    if nc == 'rRNA':
        wfile = os.path.join(args.output,args.species + '.' + nc + '_cmscan' + '.gff')
    else:
        wfile = os.path.join(args.output,args.species + '.' + nc + '.gff')
    with open(wfile,'w') as fw:
        fw.write('##gff-version 3\n')
        fw.write('\n'.join(sort_gff(gff_dict[nc])) + '\n')
