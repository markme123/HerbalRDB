# -*- coding: utf-8 -*-
"""
Created on Thu May 16 16:28:32 2019

@author: Administrator
"""

import os
import re
import argparse
import logging
import sys
import math
import time
import gzip

parser = argparse.ArgumentParser(description='cmscan ')
parser.add_argument('-r',dest='ref',type=str,help='genome fasta',default='./genome.fa')
parser.add_argument('-o',dest='output',type=str,help='output dir',default='./')
parser.add_argument('-c',dest='cmscan',type=str,help='cmscan',default='/home/software/infernal/bin/cmscan')
parser.add_argument('-Rfam_cm',dest='Rfam_cm',type=str,help='Rfam.cm',default='/home/soft_env/Rfam/Rfam.cm')
parser.add_argument('-Rfam_clanin',dest='Rfam_clanin',type=str,help='Rfam.clanin',default='/home/soft_env/Rfam/Rfam.clanin')
parser.add_argument('-s',dest='species',type=str,help='species name',default='anno')
#parser.add_argument('-size',dest='size',type=str,help='genome size',default='no')
parser.add_argument('-cm_tpye',dest='cm_tpye',type=str,help='cut_ga,cut_nc,cut_tc',choices=['cut_ga','cut_nc','cut_tc'],default='cut_ga')
parser.add_argument('-cpu',dest='cpu',type=str,help='cmscan cpus',default='6')
args = parser.parse_args()

def show_info(text):
    now_time = time.time()
    logging.info(text)
    return now_time
def run_cmd(cmd):
    logging.info(cmd)
    flag = os.system(cmd)
    if flag != 0:
        logging.error("Command fail: " + cmd)
        exit(2)
    return 0
def run_command(cmd):
    (stat1,stat2) = commands.getstatusoutput(cmd)
    return stat1,stat2
def run_time(start_time):
    spend_time = time.time() - start_time
    logging.info("Total  spend time : " + fmt_time(spend_time))
    return 0
def fmt_time(spend_time):
    spend_time = int(spend_time)
    day = 24 * 60 * 60
    hour = 60 * 60
    min = 60
    if spend_time < 60:
        return "%ds" % math.ceil(spend_time)
    elif spend_time > day:
        days = divmod(spend_time, day)
        return "%dd%s" % (int(days[0]), fmt_time(days[1]))
    elif spend_time > hour:
        hours = divmod(spend_time, hour)
        return '%dh%s' % (int(hours[0]), fmt_time(hours[1]))
    else:
        mins = divmod(spend_time, min)
        return "%dm%ds" % (int(mins[0]), math.ceil(mins[1]))
def genome_size(ref,Rfam_cm):
    genome_size = 0
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
    cm_model_num = 0
    with open(Rfam_cm,'r') as f:
        for line in f:
            if line.startswith('NAME'):
                cm_model_num += 1
    cm_model_num = cm_model_num/2
    if genome_size == 0:
        print('Please input right genome file !!!!')
        exit(1)
    else:
        #Z = total * 2 * CMmumber/106
        return format(float(genome_size*2*cm_model_num)/float(1000000),'.2f')

#if args.size != 'no' and int(args.size) != 0:
#    genome_size = args.size
#else:
#    genome_size = genome_size(args.ref,args.Rfam_cm)
logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s [line:%(lineno)d][%(levelname)s:] %(message)s',
                        datefmt='%Y-%m-%d  %H:%M:%S',
                        filename=args.output + "/cmscan.log",
                        filemode='w')
run_start = show_info("======================Run cmscan is start!======================")
#cmd = 'cd %s && %s -Z %s --%s --rfam --nohmmonly --cpu %s --tblout %s.tblout --fmt 2 --clanin %s %s %s'%(
#    args.output,args.cmscan,genome_size,args.cm_tpye,args.cpu,args.species,args.Rfam_clanin,args.Rfam_cm,args.ref)
cmd = 'cd %s && %s --%s --rfam --nohmmonly --cpu %s --tblout %s.tblout --fmt 2 --clanin %s %s %s'%(
    args.output,args.cmscan,args.cm_tpye,args.cpu,args.species,args.Rfam_clanin,args.Rfam_cm,args.ref)
print(cmd)
run_cmd(cmd)
run_time(run_start)

