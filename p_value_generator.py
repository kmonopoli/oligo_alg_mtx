#!/usr/bin/python
import csv_io
import duplex
import random_seq_generator
import parse_duplexes
import numpy as np
import math
import subprocess
import os
#input_data_file = 'normalized_nonfm_luciferase_data_cleaned_TS_101317.csv'
#input_data_file = 'nonfm_luciferase_data_cleaned_TS_101317.csv'
input_data_file = 'nonfm_luciferase_data_cleaned_TS_101317_no_GAPDH_FLUC.csv'
#input_data_file = 'normalized_nonfm_luciferase_data_cleaned_TS_101317_no_GAPDH_FLUC.csv'
#input_data_file = 'output_data_fm_cleaned-up-normalized.csv'
#input_data_file = 'output_data_luciferase_nonfm_less_than_80_removed.csv'
x = parse_duplexes.Duplex_Analysis(input_data_file)
print "WARNING: ",input_data_file," Expression %'s must be sorted lowest Expression % to highest!"

#ls = [n.twntyMer for n in x.luc_Nfm]
ls = [n.gene_region for n in x.luc_Nfm]
ls = x.get_freqsPerPos(ls)
ls = x.normalize_seq_data(ls)

# top-40 LOWEST expression % --> hits
# bottom-150 HIGHEST expression % --> fails
# nt ordering:  AUCG (same as in the weight matrix)

#temp = x.bin_by_number(x.luc_Nfm, 40, 150) ## hits, fails
#temp = x.bin_by_number(x.luc_Nfm, 90, 150) ## hits, fails
#temp = x.bin_by_number(x.luc_Nfm, 139, 150) ## hits, fails
#temp = x.bin_by_number(x.luc_Nfm,40, 88) ## hits, fails
#temp = x.bin_by_number(x.luc_Nfm, 90, 88) ## hits, fails

## for calculation without GAPDH or fLUC
#temp = x.bin_by_number(x.luc_Nfm, 33, 106) ## hits, fails 
#temp = x.bin_by_number(x.luc_Nfm, 57, 106) ## hits, fails
#temp = x.bin_by_number(x.luc_Nfm, 82, 106) ## hits, fails

#top = x.get_freqsPerPos([n.twntyMer for n in temp[0]]) ## hits
#bot = x.get_freqsPerPos([n.twntyMer for n in temp[1]]) ## fails
#top = x.get_freqsPerPos([n.eightyMer for n in temp[0]]) ## hits
#bot = x.get_freqsPerPos([n.eightyMer for n in temp[1]]) ## fails
top = x.get_freqsPerPos([n.gene_region for n in temp[0]]) ## hits
bot = x.get_freqsPerPos([n.gene_region for n in temp[1]]) ## fails

test_stats = []
i = 0
while i<len(ls[0]):
	## simulations
	P = [ls[0][i],ls[1][i],ls[2][i],ls[3][i]]
	if sum(P)!=1:
		num =  (1.0-sum(P))/4.0
		P = [p+num for p in P]

	#r_top = random_seq_generator.random_seq_gen(1,P,40)
	#r_bot= random_seq_generator.random_seq_gen(1,P,150)

	#r_top = random_seq_generator.random_seq_gen(1,P,90)
	#r_bot= random_seq_generator.random_seq_gen(1,P,150)

	#r_top = random_seq_generator.random_seq_gen(1,P,139)
	#r_bot= random_seq_generator.random_seq_gen(1,P,150)

	#r_top = random_seq_generator.random_seq_gen(1,P,90)
	#r_bot= random_seq_generator.random_seq_gen(1,P,88)

	#r_top = random_seq_generator.random_seq_gen(1,P,40)
	#r_bot= random_seq_generator.random_seq_gen(1,P,88)

	## for calculation without GAPDH or fLUC
	#r_top = random_seq_generator.random_seq_gen(1,P,33)
	#r_bot= random_seq_generator.random_seq_gen(1,P,106)

	#r_top = random_seq_generator.random_seq_gen(1,P,57)
	#r_bot= random_seq_generator.random_seq_gen(1,P,106)

	#r_top = random_seq_generator.random_seq_gen(1,P,82)
	#r_bot= random_seq_generator.random_seq_gen(1,P,106)

	## from experimental data
	t = [top[0][i],top[1][i],top[2][i],top[3][i]] 
	b = [bot[0][i],bot[1][i],bot[2][i],bot[3][i]]

	j = 0
	temp_stats = []
	while j < 4:
		## test statistic calculation

		#n1 = (float(b[j])/float(150))-(float(t[j])/float(40))
		#n2 = ((r_bot[j][2])**2/float(150)**2)/float(150)
		#n3 = ((r_top[j][2])**2/float(40)**2)/float(40)

		#n1 = (float(b[j])/float(150))-(float(t[j])/float(90))
		#n2 = ((r_bot[j][2])**2/float(150)**2)/float(150)
		#n3 = ((r_top[j][2])**2/float(90)**2)/float(90)

		#n1 = (float(b[j])/float(150))-(float(t[j])/float(139))
		#n2 = ((r_bot[j][2])**2/float(150)**2)/float(150)
		#n3 = ((r_top[j][2])**2/float(139)**2)/float(139)


		#n1 = (float(b[j])/float(88))-(float(t[j])/float(40))
		#n2 = ((r_bot[j][2])**2/float(88)**2)/float(88)
		#n3 = ((r_top[j][2])**2/float(40)**2)/float(40)

		#n1 = (float(b[j])/float(88))-(float(t[j])/float(90))
		#n2 = ((r_bot[j][2])**2/float(88)**2)/float(88)
		#n3 = ((r_top[j][2])**2/float(90)**2)/float(90)


		## for calculation without GAPDH or fLUC
		#n1 = (float(b[j])/float(106))-(float(t[j])/float(33))
		#n2 = ((r_bot[j][2])**2/float(106)**2)/float(106)
		#n3 = ((r_top[j][2])**2/float(33)**2)/float(33)

		#n1 = (float(b[j])/float(106))-(float(t[j])/float(57))
		#n2 = ((r_bot[j][2])**2/float(106)**2)/float(106)
		#n3 = ((r_top[j][2])**2/float(57)**2)/float(57)

		#n1 = (float(b[j])/float(106))-(float(t[j])/float(82))
		#n2 = ((r_bot[j][2])**2/float(106)**2)/float(106)
		#n3 = ((r_top[j][2])**2/float(82)**2)/float(82)



		n4 = n1/math.sqrt(n2+n3)
		temp_stats.append(n4)
		j+=1
	test_stats.append(temp_stats)
	i+=1
w = csv_io.CSV_Writer('test_stats.csv',['A','U','C','G'],test_stats,True)
w.write_to_file()
cwd = os.getcwd()
subprocess.call (cwd+"/p_values.R")
print "P values have been saved to pvalue.csv"
print "\a\a\a\a\a"
