import sys
import os
import glob
use_list = []
with open(sys.argv[1],'r') as f:
	for line in f:
		l = line.strip()
		if l != '':
			use_list.append(l)
with open(sys.argv[2],'w') as fw:
	for i in use_list:
		infile = glob.glob('%s.rRNA_cmscan.gff'%(i))[0]
		with open(infile,'r') as f:
			for line in f:
				l = line.strip()
				if l != '' and not l.startswith('#'):
					fw.write(l+'\n')

