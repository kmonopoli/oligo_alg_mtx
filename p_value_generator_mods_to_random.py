#!/usr/bin/python
import csv_io
import duplex
import random_seq_generator
import parse_duplexes
import numpy as np
import math
import subprocess
import os

input_data_file = 'output_data_luciferase_nonfm_less_than_80_and_DGCR8_1021_removed.csv'
x = parse_duplexes.Duplex_Analysis(input_data_file)
print "\n\n\n\n"
print "WARNING: ",input_data_file," Expression %'s must be sorted lowest Expression % to highest!"
print "\n\n\n\n"
print "WARNING: be sure P value R script has correct sample size!"
print "\n\n\n\n"
#ls = [n.sense_fm_seq for n in x.luc_Nfm]
ls = [n.antisense_fm_seq for n in x.luc_Nfm]

ls = x.get_mod_freqsPerPos(ls)
ls = x.normalize_mod_seq_data(ls)

# top-40 LOWEST expression % --> hits
# bottom-150 HIGHEST expression % --> fails
# nt ordering:  AUCG (same as in the weight matrix)

#temp = x.bin_by_number(x.luc_Nfm, 40, 150) ## hits, fails
temp = x.bin_by_number(x.luc_Nfm, 90, 230) ## hits, fails

#top = x.get_mod_freqsPerPos([n.sense_fm_seq for n in temp[0]]) ## hits
#bot = x.get_mod_freqsPerPos([n.sense_fm_seq for n in temp[1]]) ## fails
top = x.get_mod_freqsPerPos([n.antisense_fm_seq for n in temp[0]]) ## hits
bot = x.get_mod_freqsPerPos([n.antisense_fm_seq for n in temp[1]]) ## fails

test_stats = []
i = 0
while i<len(ls[0]):
	## simulations
	P = [ls[0][i],ls[1][i],ls[2][i]]#,ls[3][i]]
	rFract = 1.0/3.0
	Prand = [rFract,rFract,rFract]

	r_top = random_seq_generator.random_mod_seq_gen(1,P,90)
	r_bot= random_seq_generator.random_mod_seq_gen(1,Prand,230)

	#r_top = random_seq_generator.random_mod_seq_gen(1,Prand,90)
	#r_bot= random_seq_generator.random_mod_seq_gen(1,P,230)

	## from experimental data
	t = [top[0][i],top[1][i],top[2][i]]#,top[3][i]] 
	b = [bot[0][i],bot[1][i],bot[2][i]]#,bot[3][i]]

	j = 0
	temp_stats = []
	while j < 3:#4:
		## test statistic calculation
		n1 = (float(b[j])/float(230))-(float(t[j])/float(90))
		n2 = ((r_bot[j][2])**2/float(230)**2)/float(230)
		n3 = ((r_top[j][2])**2/float(90)**2)/float(90)


		#n1 = (float(b[j])/float(150))-(float(t[j])/float(40))
		#n2 = ((r_bot[j][2])**2/float(150)**2)/float(150)
		#n3 = ((r_top[j][2])**2/float(40)**2)/float(40)

		if (n2+n3) == 0: # top prevent divide by 0 error
			n4 = 0
		else:
			n4 = n1/math.sqrt(n2+n3)
		temp_stats.append(n4)
		j+=1
	test_stats.append(temp_stats)
	i+=1
w = csv_io.CSV_Writer('test_stats_mod.csv',['f','m','non'],test_stats,True)
w.write_to_file()
cwd = os.getcwd()
subprocess.call (cwd+"/p_values_mod.R")
print "P values have been saved to pvalue_mod.csv"
print "\a\a\a\a\a"
