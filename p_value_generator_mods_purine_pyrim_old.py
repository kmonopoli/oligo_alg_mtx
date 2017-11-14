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
print "\n\n"
print "WARNING: ",input_data_file," Expression %'s must be sorted lowest Expression % to highest!"
print "\n\n"
print "WARNING: be sure P value R script has correct sample size!"
print "\n\n"
#ls = [n.sense_fm_seq for n in x.luc_Nfm]
ls_1 = [n.twntyMer for n in x.luc_Nfm]
ls_1 = x.make_20mer_antisense(ls_1)
ls_1 = x.get_freqsPerPos(ls_1)
ls_1 = x.normalize_seq_data(ls_1)

ls_2 = [n.antisense_fm_seq for n in x.luc_Nfm]
ls_2 = x.get_mod_freqsPerPos_purines_vs_pyrimidines(ls_2)
ls_2 = x.normalize_mod_seq_data_purine_pyrimidine(ls_2)


# top-40 LOWEST expression % --> hits
# bottom-150 HIGHEST expression % --> fails
# nt ordering:  AUCG (same as in the weight matrix)

temp = x.bin_by_number(x.luc_Nfm, 40, 150) ## hits, fails

top_1 = x.get_freqsPerPos(x.make_20mer_antisense([n.twntyMer for n in temp[0]])) ## hits
bot_1 = x.get_freqsPerPos(x.make_20mer_antisense([n.twntyMer for n in temp[1]])) ## fails

top_2 = x.get_mod_freqsPerPos_purines_vs_pyrimidines([n.antisense_fm_seq for n in temp[0]]) ## hits
bot_2 = x.get_mod_freqsPerPos_purines_vs_pyrimidines([n.antisense_fm_seq for n in temp[1]]) ## fails

# need to calculate two P's for modifications at each position depending on if it is a purine or pyrimidine (or else probabilities won't sum to 1 and random sequence generator won't work)

test_stats = []
i = 0
while i<len(ls_1[0]):
	## simulations
	P_1= [ls_2[0][0][i],ls_2[0][1][i],ls_2[0][2][i]] # pyrimidines CU (f,m,non)
	P_2= [ls_2[1][0][i],ls_2[1][1][i],ls_2[1][2][i]] # purines AG (f,m,non)
	P_3 = [ls_1[0][i],ls_1[1][i],ls_1[2][i],ls_1[3][i]] # sequences (A,U,C,G)

	if sum(P_2) == 0:# case where probabilities of modifications for purines = 0, set f to 1 
		if i == 0:
			if P_3[0] == 0 and P_3[3] == 0:	#checks that P of A or G at this position is 0 to ensure this has no effect on results
				P_2 = [1.0,0.0,0.0]
			else:
				print "WARNING: probabilities of modifications of purines is zero while probabilities of purines is not!"
		else:
			print "WARNING: probabilities of modifications of purines is zero at position ",i

	if (sum(P_1) < 0.99) or (sum(P_2) < 0.99) or (sum(P_3) < 0.99):
		print "WARNING: probabilities do not sum to one for at least one list"

	r_top = random_seq_generator.random_seq_and_mod_gen(1,P_3,P_1,P_2,40)
	r_bot = random_seq_generator.random_seq_and_mod_gen(1,P_3,P_1,P_2,150)


	## from experimental data
	t =[top_1[0][i], top_1[1][i], top_1[2][i], top_1[3][i],
		top_2[0][0][i], top_2[0][1][i], top_2[0][2][i],
		top_2[1][0][i], top_2[1][1][i], top_2[1][2][i]]

	b =[bot_1[0][i], bot_1[1][i], bot_1[2][i], bot_1[3][i],
		bot_2[0][0][i], bot_2[0][1][i], bot_2[0][2][i],
		bot_2[1][0][i], bot_2[1][1][i], bot_2[1][2][i]]

	j = 0
	temp_stats = []
	while j < 10:#3:#4:
		## test statistic calculation
		n1 = (float(b[j])/float(150))-(float(t[j])/float(40))
		#n1 = (float(b[j])/float(150))-(float(t[j])/float(40))
		n2 = ((r_bot[j][2])**2/float(150)**2)/float(150)
		#n2 = ((r_bot[j][2])**2/float(150)**2)/float(150)
		n3 = ((r_top[j][2])**2/float(40)**2)/float(40)
		if (n2+n3) == 0: # top prevent divide by 0 error
			n4 = 0
		else:
			n4 = n1/math.sqrt(n2+n3)
		temp_stats.append(n4)
		j+=1
	test_stats.append(temp_stats)
	i+=1
w = csv_io.CSV_Writer('test_stats_mod_and_seq.csv',['A','U','C','G','pyr_f','pyr_m','pyr_non','pur_f','pur_m','pur_non'],test_stats,True)
w.write_to_file()
cwd = os.getcwd()
subprocess.call (cwd+"/p_values_mod_and_seq.R")
print "P values have been saved to pvalue_mod_and_seq.csv"
print "\a\a\a\a\a"
