#!/usr/bin/python
import csv_io
import duplex
import random_seq_generator
import parse_duplexes
import numpy as np
import math
import os


# read in fluch probabilities
probabilities_fluch_file = 'probabilities_FLUCH_for_random_p_value_calc_old.csv'
fluch_prob_reader= csv_io.CSV_Reader(probabilities_fluch_file)
fluch_prob_reader.read_file()
Pfluch = fluch_prob_reader.file_data
Pfluch = [[float(n.replace('\n','')) for n in q] for q in Pfluch]

#seqs_t = []
#seqs_b = []
seqs = []
i = 0
while i<45:#len(Pfluch):
	## simulations
	#Prand = [0.265035596243,0.250894084916,0.240274032276,0.243796286564]
	Prand = Pfluch[i]
	if sum(Prand)!=1:
		if (sum(Prand) < 0.99) or (sum(Prand) > 1.01):
			print "WARNING: FLUCH probabilities at position ",i," do not sum to 1 !"
			print "sum to: ",sum(Prand)
		num =  (1.0-sum(Prand))/4.0
		Prand = [p+num for p in Prand]
	print Prand

	rand = random_seq_generator.make_rand_seqs(1,Prand,216)
#	rand_t = random_seq_generator.random_seq_gen(1,Prand,57)
#	rand_t = [x[0]/57.0 for x in rand_t]
#	seqs_t.append(rand_t)

#	rand_b = random_seq_generator.random_seq_gen(1,Prand,106)
#	rand_b = [x[0]/106.0 for x in rand_b]
#	seqs_b.append(rand_b)
	seqs.append(rand)
	i+=1
#seqs_t.append([0.0,0.0,0.0,0.0])
#seqs = seqs_t + seqs_b

#w = csv_io.CSV_Writer('rand_seqs_refseq.csv',["A","U","C","G"],seqs,True)
w = csv_io.CSV_Writer('rand_seqs_refseq.csv',list(range(1,len(seqs)+1)),seqs,False)
#w = csv_io.CSV_Writer('rand_seqs_refseq.csv',list(range(1,45+1)),seqs,True)
w.write_to_file()
