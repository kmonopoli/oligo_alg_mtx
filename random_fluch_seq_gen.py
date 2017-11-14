###### 	OLD ####### USE random_seq_gen_from_data_freqs.py ###

#!/usr/bin/python
import csv_io
import duplex
import random_seq_generator
import parse_duplexes
import numpy as np
import math
import os


def get_freqsPerPos(ls):
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


def normalize_seq_data(ls):
	s= ls[0][0]+ls[1][0]+ls[2][0]+ls[3][0]
	return [[float(n)/s for n in k] for k in ls]


# read in fluch probabilities
probabilities_fluch_file = 'output/probabilities_FLUCH_for_random_p_value_calc_110617.csv'
fluch_prob_reader= csv_io.CSV_Reader(probabilities_fluch_file)
fluch_prob_reader.read_file()
Pfluch = fluch_prob_reader.file_data
#Pfluch = [[float(n.replace('\n','')) for n in q] for q in Pfluch]
Pfluch = [q[:4] for q in Pfluch]
Pfluch = [[float(n) for n in q] for q in Pfluch]
dataLs = []
j = 0
while j<25:
	seqs = []
	i = 0
	while i<45:#len(Pfluch):
		## simulations
		Prand = Pfluch[i]
		if sum(Prand)!=1:
			if (sum(Prand) < 0.99) or (sum(Prand) > 1.01):
				print "WARNING: FLUCH probabilities at position ",i," do not sum to 1 !"
				print "sum to: ",sum(Prand)
			num =  (1.0-sum(Prand))/4.0
			Prand = [p+num for p in Prand]

		#rand = random_seq_generator.make_rand_seqs(1,Prand,216)
		rand = random_seq_generator.make_rand_seqs(1,Prand,158)
		seqs.append(rand)
		i+=1
	#split into top and bottom bins
	#top = 57
	#bottom = 106
	#seqs = [[x[q] for x in seqs] for q in list(range(0,216))] # transposes data 
	seqs = [[x[q] for x in seqs] for q in list(range(0,158))] # transposes data 
	#top = seqs[:57]
	top = seqs[:51]
	#bottom = seqs[(216-106):]
	bottom = seqs[(158-57):]
	top = normalize_seq_data(get_freqsPerPos(top))
	bottom = normalize_seq_data(get_freqsPerPos(bottom))
	diff_A = [ top[0][i]-bottom[0][i] for i in list(range(0,45))]
	diff_U = [ top[1][i]-bottom[1][i] for i in list(range(0,45))]
	diff_C = [ top[2][i]-bottom[2][i] for i in list(range(0,45))]
	diff_G = [ top[3][i]-bottom[3][i] for i in list(range(0,45))]
	ls = [diff_A,diff_U,diff_C,diff_G,[""]*45]

	dataLs+=ls
	j+=1
	



w = csv_io.CSV_Writer('rand_seqs_refseq.csv',list(range(1,45+1)),dataLs,True)
w.write_to_file()
