#!/usr/bin/python
import csv_io
import duplex
import random_seq_generator
import parse_duplexes
import numpy as np
import math
import subprocess
import os

input_data_file = 'csv_dbcomp/normalized_nonfm_luciferase_data_cleaned_TS_101317.csv'#'output_data_diana.csv'
x = parse_duplexes.Duplex_Analysis(input_data_file)
print "WARNING: ",input_data_file," Expression %'s must be sorted lowest Expression % to highest!"

ls = [n.sense_fm_seq for n in x.luc_Nfm]
#ls = [n.antisense_fm_seq for n in x.luc_Nfm]
wts = [q[1:4] for q in csv_io.CSV_Reader('output/weight_matrixes_for_paper/mod_weight_matrix_top90_bot150_luc_nonfm_sense.csv').file_data]
wts = [[int(q1) for q1 in q2] for q2 in wts] #[f,m,non]
def score_seq(seq,weights):
	scr = 0
	if len(seq) != len(weights):
		print "ERROR: sequence length is ",len(seq)," but weights length is ",len(weights)
		quit()
	i = 0
	while i<len(seq):
		s = seq[i].replace('G','').replace('A','').replace('U','').replace('C','').replace('1','').replace('P','')
		if(s == "f"):	
			scr += weights[i][0]
		elif(s == "m"):
			scr += weights[i][1]
		elif(s == ""):
			scr += weights[i][2]
		else:
			print "WARNING: value ",s," is not a known modification so scoring is possibly inaccurate"
		i+=1
	return scr

scores = []
for se in ls:
	scores.append(score_seq(se,wts))

id_ls = [q3.iD for q3 in x.luc_Nfm]

data_ls = []
i = 0
while i<len(id_ls):
	data_ls.append([id_ls[i],scores[i]])
	i+=1

w = csv_io.CSV_Writer('mod_scores.csv',['iD','modScore'],data_ls,True)
w.write_to_file()
#ls = x.get_mod_freqsPerPos(ls)
#ls = x.normalize_mod_seq_data(ls)
