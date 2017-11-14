#!/usr/bin/python
import csv_io
import duplex
import random_seq_generator
import parse_duplexes
import numpy as np
import math
import subprocess
import os

#input_data_file = 'output_data_fm_cleaned-up-normalized.csv'
input_data_file = 'output_data_fm_unclean-normalized_for_scoring.csv'

x = parse_duplexes.Duplex_Analysis(input_data_file)

#ls = [n.sense_fm_seq for n in x.luc_Nfm]
#ls = [n.antisense_fm_seq for n in x.luc_Nfm]

ls = [n.gene_region for n in x.fm if len(n.gene_region) == 45]
#ls = [n.gene_region for n in x.all_duplexes if len(n.gene_region) == 45]
#ls = [n.eightyMer for n in x.fm if len(n.eightyMer) ==80]
print len(ls)

#dir1 = 'fm/top40_bottom70/weight_matrix_fm_top40_bot70.csv'
#dir1 = 'weight_matrix_FLUCH.csv'
dir1 = 'weight_matrix_top40_bot150_AVG10_10.csv'

wts = [q[1:5] for q in csv_io.CSV_Reader('output/'+dir1).file_data]
wts = [[int(q1) for q1 in q2] for q2 in wts] #[a,u,c,g]
def score_seq(seq,weights):
	scr = 0
	if len(seq) != len(weights):
		print "ERROR: sequence length is ",len(seq)," but weights length is ",len(weights)
		quit()
	i = 0
	while i<len(seq):
		s = seq[i]
		if(s == "A"):	
			scr += weights[i][0]
		elif(s == "U"):
			scr += weights[i][1]
		elif(s == "C"):
			scr += weights[i][2]
		elif(s == "G"):
			scr += weights[i][3]
		else:
			print "WARNING: value ",s," is not a known modification so scoring is possibly inaccurate"
		i+=1
	return scr


id_ls = [q3.iD for q3 in x.fm if len(q3.gene_region) == 45]
#id_ls = [q3.iD for q3 in x.all_duplexes if len(q3.gene_region) == 45]
#id_ls = [q3.iD for q3 in x.fm if len(q3.eightyMer) == 80]


scores = []
for se in ls:
	scores.append(score_seq(se,wts))

data_ls = []
i = 0
while i<len(id_ls):
	data_ls.append([id_ls[i],scores[i]])
	i+=1

w = csv_io.CSV_Writer('seq_scores_fm.csv',['iD','Score'],data_ls,True)
w.write_to_file()
#ls = x.get_mod_freqsPerPos(ls)
#ls = x.normalize_mod_seq_data(ls)
