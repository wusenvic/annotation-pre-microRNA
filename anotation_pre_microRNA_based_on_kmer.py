"""
first to read the mirtable generated from vsearch with parameters as below:
vsearch --usearch_global miR/SRR5804875_1.fa --db ../reference/hairpin_hsa.fa --id 0.97 --otutabout otu_table3.txt --matched match3.fa --blast6out tnbc_miR.blast3 --mincols 15 --matched SRR10248408_b20.fa

first: mirtable form vsearch and get transted
second: blast result from vsearch
third: new table of the mirtable
fourth: sequense file selected based on the blast table
fifth: mature miR index
sixth: premature miR index
"""

def generate_kmers(sequence, k):
	kmers = []
	for i in range(len(sequence) - k + 1):
		kmers.append(sequence[i:i+k])
	return kmers

import os
import sys
import re

raw_table = {}
raw_seq_number = 0
b = {}
cal_factors = {}
new_table = {}
header = []
raw_seq = {}
M_index = {}
pre_index = {}

final_index = {}

query_seq = {}

match_seq = {}
file_size_line = {}
mir_array = {}



with open(sys.argv[5],'r') as index_file:
	for line in index_file:
		if '>' in line:
			if "-3p" in line or "-5p" in line:
				name = line.split(' ')[0][1:-3]
			else:
				name = line.split(' ')[0][1:]
			if name in M_index.keys():
				next
			else:
				M_index[name] = 0
		else:
			leng_M = len(line.strip())
			if leng_M >= M_index[name]:
				M_index[name]=leng_M
			else:
				next

with open(sys.argv[6],'r') as index_file:
	for line in index_file:
		if line[0] == '>':
			name = line.split(' ')[0][1:]
			pre_index[name] = [""]
			pre_index[name][0] = ""
			mir_array[name] = 0
		else:
			pre_index[name][0] += line.strip()
	for pre_name, pre_seqs in pre_index.items():
		if pre_name not in M_index.keys():
#			print(pre_name)
			new_name_m = '-'.join(pre_name.split('-')[:-1])
			kmer_seq = generate_kmers(pre_seqs[0], M_index[new_name_m])
		else:
			kmer_seq = generate_kmers(pre_seqs[0], M_index[pre_name])
		pre_seqs.extend(kmer_seq)

with open(sys.argv[1],'r') as mirT:
	for line in mirT:
		if line[0] != '#':
			raw_table[line.split(' ')[0]] = line.strip().split(' ')[1:]
			for i in line.strip().split(' ')[1:]:
				raw_seq_number += int(i)
		else:
			header = line.strip().split(' ')[1:]


with open(sys.argv[4]) as seqs_to_filter:
	for lines in seqs_to_filter:
		lines.strip()
		if lines[0] == '>':
			seq_name  = lines[1:].strip()
#			print(seq_name)
			raw_seq[seq_name] = ""
		else:
			raw_seq[seq_name] += lines.strip()
with open(sys.argv[2],'r') as blastT:
	for line in blastT:
		l = line.strip().split('\t')
		sample_seq_name = l[0]
		sample_name = sample_seq_name.split('.')[0]
		mir_index = header.index(l[1])
		if l[1] not in M_index.keys():
			new_char = '-'.join(pre_name.split('-')[:-1])
			raw_seq_k = generate_kmers(raw_seq[sample_seq_name],M_index[new_char])
		else:
			raw_seq_k = generate_kmers(raw_seq[sample_seq_name],M_index[l[1]])
		num = 0
		k_index = []
		for kmer in raw_seq_k:
			if kmer in pre_index[l[1]][1:]:
				num += 1
				k_index.append(raw_seq_k.index(kmer))
		b[sample_seq_name] = [sample_seq_name,sample_name,mir_index,num,len(pre_index[l[1]][1:]),k_index]
new_table_count = raw_seq_number

critia =3
new_table = raw_table
new_table_count = raw_seq_number
full_k = len(pre_index[l[1]][1:])
for sample_seq_name, vals in b.items():
	if vals[3] < critia:
		new_table[vals[1]][vals[2]] = str(int(new_table[vals[1]][vals[2]])-1)
		new_table_count = new_table_count - 1
	elif vals[3] > vals[4]:
		new_table[vals[1]][vals[2]] = str(int(new_table[vals[1]][vals[2]])+1)
	else:
		if vals[-1][-1] - vals[-1][-2] > vals[3] -1 :
			new_table[vals[1]][vals[2]] = str(int(new_table[vals[1]][vals[2]])-1)
		else:
			new_table[vals[1]][vals[2]] = str(int(new_table[vals[1]][vals[2]]))
with open(sys.argv[3],'w') as new:
	new.write('samples\t'+'\t'.join(header)+'\n')
	for samples in new_table.keys():
		new.write(samples+'\t'+'\t'.join(new_table[samples])+'\n')

