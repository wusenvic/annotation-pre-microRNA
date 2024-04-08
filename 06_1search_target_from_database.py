#! usr/bin/env python3
import os
import sys
import re
import requests
url = 'http://www.mirdb.org/cgi-bin/search.cgi?searchType=miRNA&searchBox='
target_scan_url = 'https://www.targetscan.org/cgi-bin/targetscan/vert_80/targetscan.cgi?species=Human&gid=&mir_sc=&mir_c=&mir_nc=&mir_vnc=&mirg='
suffix = '&full=1'
pattern = re.compile(r'>(\w+| \w+)<')
target_scan_pattern = re.compile('>(\S+)<')
target = {}
f = open(sys.argv[3],'w')
with open(sys.argv[1],'r') as mir_list:
	for lines in mir_list:
		print('processing microRNA:\t'+lines)
		target[lines.strip()] = set()
		ful_url = url+lines.strip()+suffix
		response = str(requests.get(ful_url).content).split('<tr')
		for i in response:
			b = pattern.findall(i)
			if len(b) > 3 and int(b[2]) >= 90:
#				print(b[2])
				f.write(lines.strip()+'\t'+'\t'.join(b)+'\n')
				target[lines.strip()].add(b[-1])
		ful_targetscan_url = target_scan_url+lines.strip()[4:]
		tsc_response3 = str(requests.get(ful_targetscan_url).content)
		if 'You' not in tsc_response3:
			tsc_response =tsc_response3.split('\\n')
			c = []
			for i in tsc_response:
				b = target_scan_pattern.findall(i)
				if len(b) > 5:
					n = ''.join(b)
					n = re.sub(r'</td>','\t',n)
					x = re.sub(r'<\S+?>','',n)
					x = re.sub(r'\t\t','\t',x)
					if x.count('\t')>3:
						c.append(x)
			index = 0
			print(len(c))
			l = c[0].split('\t')
			for i in l:
				if re.search('^-',i):
					index = l.index(i)
	#				print(i)
			for info in c:
				if info[:5]=='total':
					info = re.sub(r'total8mer7mer-m87mer-','',info)
				f.write(info+'\n')
				t = info.split('\t')
	#			print(t)
				if float(t[index]) <= -0.50:
					target[lines.strip()].add(t[0])
f.close()
with open(sys.argv[2],'w') as target_file:
	for mir, target_genes in target.items():
		c = '\t'.join(list(target_genes))
		target_file.write(mir+'\t'+c+'\n')
