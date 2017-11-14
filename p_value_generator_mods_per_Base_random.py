#!/usr/bin/python
import csv_io
import duplex
import random_seq_generator
import parse_duplexes
import numpy as np
import math
import subprocess
import os
input_data_file = 'csv_dbcomp/normalized_nonfm_luciferase_data_cleaned_TS_101317.csv'
x = parse_duplexes.Duplex_Analysis(input_data_file)
print "WARNING: ",input_data_file," Expression %'s must be sorted lowest Expression % to highest!"
ls = [n.antisense_fm_seq for n in x.luc_Nfm]
ls = x.get_mod_freqsPerPos_per_base(ls)
ls = x.normalize_mod_seq_data_per_base(ls)

# read in fluch probabilities
probabilities_fluch_file = 'mod_per_base_probabilities_as.csv'
fluch_prob_reader= csv_io.CSV_Reader(probabilities_fluch_file)
fluch_prob_reader.read_file()
Pfluch = fluch_prob_reader.file_data
Pfluch = [l[1:13] for l in Pfluch]
Pfluch = [[float(n.replace('\n','')) for n in q] for q in Pfluch]

# top-40 LOWEST expression % --> hits
# bottom-150 HIGHEST expression % --> fails
# nt ordering:  AUCG (same as in the weight matrix)

#temp = x.bin_by_number(x.luc_Nfm, 40, 150) ## hits, fails
temp = x.bin_by_number(x.luc_Nfm, 90, 150) ## hits, fails
#temp = x.bin_by_number(x.luc_Nfm, 139, 150) ## hits, fails

top_2 = x.get_mod_freqsPerPos_per_base([n.antisense_fm_seq for n in temp[0]]) ## hits
bot_2 = x.get_mod_freqsPerPos_per_base([n.antisense_fm_seq for n in temp[1]]) ## fails

test_stats = []
i = 0
while i<len(ls[0]):
	## simulations

	P= [
		ls[0][i],ls[1][i],ls[2][i]
		,ls[3][i],ls[4][i],ls[5][i]
		,ls[6][i],ls[7][i],ls[8][i]
		,ls[9][i],ls[10][i],ls[11][i]
		]
	if (sum(P) != 1):
		print "WARNING: does not sum to 1: ",sum(P)
		n = (1.0-sum(P))
		k =0
		while(P[k]==0):
			k+=1
		P[k]+=n

	Prand = Pfluch[i] 
	if sum(Prand)!=1:
		n = (1.0-sum(Prand))
		k =0
		while(Prand[k]==0):
			k+=1
		Prand[k]+=n


	r_top = random_seq_generator.random_mod_gen_per_base(1,Prand,90)
	r_bot= random_seq_generator.random_mod_gen_per_base(1,Prand,150)




	## from experimental data
	t =[
		top_2[0][i], top_2[1][i], top_2[2][i],
		top_2[3][i], top_2[4][i], top_2[5][i],
		top_2[6][i], top_2[7][i], top_2[8][i],
		top_2[9][i], top_2[10][i], top_2[11][i]]
	b =[
		bot_2[0][i], bot_2[1][i], bot_2[2][i],
		bot_2[3][i], bot_2[4][i], bot_2[5][i],
		bot_2[6][i], bot_2[7][i], bot_2[8][i],
		bot_2[9][i], bot_2[10][i], bot_2[11][i]]

	
	j = 0
	temp_stats = []
	while j < 12:
		## test statistic calculation

		#n1 = (float(b[j])/float(150))-(float(t[j])/float(40))
		#n2 = ((r_bot[j][2])**2/float(150)**2)/float(150)
		#n3 = ((r_top[j][2])**2/float(40)**2)/float(40)

		n1 = (float(b[j])/float(150))-(float(t[j])/float(90))
		n2 = ((r_bot[j][2])**2/float(150)**2)/float(150)
		n3 = ((r_top[j][2])**2/float(90)**2)/float(90)

		#n1 = (float(b[j])/float(150))-(float(t[j])/float(139))
		#n2 = ((r_bot[j][2])**2/float(150)**2)/float(150)
		#n3 = ((r_top[j][2])**2/float(139)**2)/float(139)

		
		# to prevent divide by 0
		if (n2+n3) == 0: 
			n4 = 0
		else:
			n4 = n1/math.sqrt(n2+n3)

		temp_stats.append(n4)
		j+=1
	test_stats.append(temp_stats)
	i+=1
#w = csv_io.CSV_Writer('test_stats.csv',['f','m','non'],test_stats,True)

w = csv_io.CSV_Writer('test_stats_mod_per_base.csv',['a_f','a_m','a_non','u_f','u_m','u_non','c_f','c_m','c_non','g_f','g_m','g_non',],test_stats,True)
w.write_to_file()
cwd = os.getcwd()
subprocess.call (cwd+"/p_values_mod_per_base.R")
print "P values have been saved to pvalue_mod_per_base.csv"
print "\a\a\a\a\a"

