#!/usr/bin/python
import csv_io
import duplex
import random_seq_generator
import parse_duplexes
import numpy as np
import math
import os


def get_freqsPerPos1(ls):
	freqA = [0]*len(ls[0])
	freqU = [0]*len(ls[0])
	freqC = [0]*len(ls[0])
	freqG = [0]*len(ls[0])
	for seq in ls:
		i = 0
		while i < len(seq):
			if seq[i] == 'A':
				freqA[i]+=1

			elif seq[i] == 'U':
				freqU[i]+=1

			elif seq[i] == 'C':
				freqC[i]+=1

			elif seq[i] == 'G':
				freqG[i]+=1
			else:
				print "ERROR: sequence ",ls ," contains non AUCG base: ",seq[i]
			i+=1
	return [freqA,freqU,freqC, freqG]


def normalize_seq_data1(ls):
	s= ls[0][0]+ls[1][0]+ls[2][0]+ls[3][0]
	return [[float(n)/s for n in k] for k in ls]



#input_data_file = 'normalized_nonfm_luciferase_data_cleaned_TS_101317.csv'
#input_data_file = 'nonfm_luciferase_data_cleaned_TS_101317.csv'
#input_data_file = 'nonfm_luciferase_data_cleaned_TS_101317_no_GAPDH_FLUC.csv'
input_data_file = 'sd_data_clean.csv'
#input_data_file = 'HTT_data_non_FLUCH_for_random_seq_analysis.csv'
#input_data_file = 'normalized_nonfm_luciferase_data_cleaned_TS_101317_no_GAPDH_FLUC.csv'
#input_data_file = 'output_data_fm_cleaned-up-normalized.csv'
#input_data_file = 'output_data_luciferase_nonfm_less_than_80_removed.csv'
#input_data_file = 'qPCR_nonfm_normalized_data.csv'

x = parse_duplexes.Duplex_Analysis(input_data_file)


#ls = [n.gene_region for n in x.luc_Nfm]
ls = [n.gene_region for n in x.qpcr_Nfm]
#ls = [n.eightyMer for n in x.qpcr_Nfm]
num_seqs = len(ls)
print num_seqs
ls = x.get_freqsPerPos(ls)
ls = x.normalize_seq_data(ls)
#bcut = 57 #bottom cutoff
#tcut = 106 #top cutoff
bcut = 57 #bottom cutoff
tcut = 51 #top cutoff
#bcut = 14 #bottom cutoff
#tcut = 49 #top cutoff
dataLs = []
j = 0
while j<num_seqs:
	seqs = []
	i = 0
	while i<45:#80:#45:#len(Pfluch):
		## simulations
		P = [ls[0][i],ls[1][i],ls[2][i],ls[3][i]]
		if sum(P)!=1:
			num =  (1.0-sum(P))/4.0
			P = [p+num for p in P]
		rand = random_seq_generator.make_rand_seqs(1,P,num_seqs)
		seqs.append(rand)
		i+=1

	#split into top and bottom bins
	seqs = [[z[q] for z in seqs] for q in list(range(0,num_seqs))] # transposes data 
	top = seqs[:bcut]
	bottom = seqs[(num_seqs-tcut):]
	top = normalize_seq_data1(get_freqsPerPos1(top))
	bottom = normalize_seq_data1(get_freqsPerPos1(bottom))
	diff_A = [ top[0][i]-bottom[0][i] for i in list(range(0,45))]
	diff_U = [ top[1][i]-bottom[1][i] for i in list(range(0,45))]
	#diff_A = [ top[0][i]-bottom[0][i] for i in list(range(0,80))]
	#diff_U = [ top[1][i]-bottom[1][i] for i in list(range(0,80))]
	## GC will just be opposite of AU so no need to calculate

	diff_AU = [sum(diff_A[i:i+5])+sum(diff_U[i:i+5]) for i in list(range(0,46-5))]
	#diff_AU = [sum(diff_A[i:i+5])+sum(diff_U[i:i+5]) for i in list(range(0,81-5))]

	dls = [diff_AU]#,diff_CG]#,[""]*45]

	dataLs+=dls
	j+=1

dataLs = [[z[q] for z in dataLs] for q in list(range(0,46-5))] # transposes data 
#dataLs = [[z[q] for z in dataLs] for q in list(range(0,81-5))] # transposes data 
for q in dataLs:
	list.sort(q)

w = csv_io.CSV_Writer('rand_seqs_refseq.csv',list(range(1,46-5)),dataLs,False)
#w = csv_io.CSV_Writer('rand_seqs_refseq.csv',list(range(1,81-5)),dataLs,False)
w.write_to_file()



