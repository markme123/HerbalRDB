import sys
import os
if not os.path.exists(sys.argv[1]):
	os.system('touch %s'%(sys.argv[2]))
	exit(0)
with open(sys.argv[1],'r') as f,open(sys.argv[2],'w') as fw:
	for line in f:
		l = line.strip()
		if l != '' and not l.startswith('#'):
				ll = l.split('\t')
				if ll[-1] == '8s_rRNA':
					ll[-1] = 'Annotation=5.8S ribosomal RNA'
					fw.write('\t'.join(ll)+'\n')
				elif ll[-1] == '28s_rRNA':
					ll[-1] = 'Annotation=28S ribosomal RNA'
					fw.write('\t'.join(ll)+'\n')
				elif ll[-1] == '18s_rRNA':
					ll[-1] = 'Annotation=18S ribosomal RNA'
					fw.write('\t'.join(ll)+'\n')
				else:
					pass
