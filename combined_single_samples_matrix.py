#! /usr/bin/env python3

import os
import sys

new_matrix = {}
new_header = []
com_matrix = {}

def read_matrix(line):
	microRNA= {}
	with open(line) as raw_file:
		for lines in raw_file:
			if lines[0:3] !='hsa':
				header = '\t'.join(lines.strip().split('\t')[1:])
			else:
				microRNA[lines.split('\t')[0]] = lines.strip().split('\t')[1:]
	return [header,microRNA]

def combine(new_matrix, single_matrix,header):
	inter_media = {}
	if len(new_matrix) == 0:
		new_matrix = single_matrix[1]
	else:
		a=['0']
		len_a = len(header)-1

		for miR in single_matrix[1].keys():
			if miR in new_matrix.keys():
				inter_media[miR] = new_matrix[miR]+single_matrix[1][miR]
			else:
				pre_num = a*len_a

				type(pre_num)
				inter_media[miR] = pre_num+single_matrix[1][miR]
		for miR in new_matrix.keys():
			if miR not in single_matrix[1].keys():
				pre_num = a*len_b

				inter_media[miR] = new_matrix[miR]+a
			else:
				next

		new_matrix = inter_media
	return new_matrix
with open(sys.argv[1],'r') as file_list:
	for line in file_list:
		first_matrix = read_matrix(line.strip())
		with open(line.strip()+'single.txt','w') as new:
			new.write('microRNA\t'+first_matrix[0]+'\n')
			for a,b in first_matrix[1].items():
				new.write(a+'\t'+'\t'.join(b)+'\n')

		len_b = len(first_matrix[0])
		print(len_b)

		new_header.append(first_matrix[0])
		new_matrix = combine(new_matrix,first_matrix,new_header)

 
print(new_header)
with open(sys.argv[2],'w') as final_file:
		final_file.write('microRNA'+'\t'+'\t'.join(new_header)+'\n')
		for miR,exp in new_matrix.items():
			line_sum = 0
			for num in exp:
				line_sum += int(num)
			if line_sum != 0:
				final_file.write(miR+'\t'+'\t'.join(exp)+'\n')
			else:
				next
