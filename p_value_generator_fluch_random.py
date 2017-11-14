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
input_data_file = 'normalized_nonfm_luciferase_data_cleaned_TS_101317_no_GAPDH_FLUC.csv'
x = parse_duplexes.Duplex_Analysis(input_data_file)
print "WARNING: ",input_data_file," Expression %'s must be sorted lowest Expression % to highest!"
ls = [n.gene_region for n in x.luc_Nfm]
ls = x.get_freqsPerPos(ls)
ls = x.normalize_seq_data(ls)

# read in fluch probabilities
probabilities_fluch_file = 'probabilities_FLUCH_for_random_p_value_calc.csv'
fluch_prob_reader= csv_io.CSV_Reader(probabilities_fluch_file)
fluch_prob_reader.read_file()
Pfluch = fluch_prob_reader.file_data
Pfluch = [[float(n.replace('\n','')) for n in q] for q in Pfluch]

# top-40 LOWEST expression % --> hits
# bottom-150 HIGHEST expression % --> fails
# nt ordering:  AUCG (same as in the weight matrix)

#temp = x.bin_by_number(x.luc_Nfm, 40, 150) ## hits, fails
#temp = x.bin_by_number(x.luc_Nfm, 90, 150) ## hits, fails
#temp = x.bin_by_number(x.luc_Nfm, 139, 150) ## hits, fails

## for calculation without GAPDH or fLUC
#temp = x.bin_by_number(x.luc_Nfm, 40, 98) ## hits, fails
#temp = x.bin_by_number(x.luc_Nfm, 58, 98) ## hits, fails
temp = x.bin_by_number(x.luc_Nfm, 86, 98) ## hits, fails

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

	Prand = Pfluch[i] 
	if sum(Prand)!=1:
		if (sum(Prand) < 0.99) or (sum(Prand) > 1.01):
			print "WARNING: FLUCH probabilities at position ",i," do not sum to 1 !"
			print "sum to: ",sum(Prand)
		num =  (1.0-sum(Prand))/4.0
		Prand = [p+num for p in Prand]

	#r_top = random_seq_generator.random_seq_gen(1,P,40)
	#r_bot= random_seq_generator.random_seq_gen(1,Prand,150)

	#r_top = random_seq_generator.random_seq_gen(1,Prand,40)
	#r_bot= random_seq_generator.random_seq_gen(1,P,150)

	#r_top = random_seq_generator.random_seq_gen(1,Prand,40)
	#r_bot= random_seq_generator.random_seq_gen(1,Prand,150)


	#r_top = random_seq_generator.random_seq_gen(1,P,90)
	#r_bot= random_seq_generator.random_seq_gen(1,Prand,150)

	#r_top = random_seq_generator.random_seq_gen(1,Prand,90)
	#r_bot= random_seq_generator.random_seq_gen(1,P,150)

	#r_top = random_seq_generator.random_seq_gen(1,Prand,90)
	#r_bot= random_seq_generator.random_seq_gen(1,Prand,150)


	#r_top = random_seq_generator.random_seq_gen(1,P,139)
	#r_bot= random_seq_generator.random_seq_gen(1,Prand,150)

	#r_top = random_seq_generator.random_seq_gen(1,Prand,139)
	#r_bot= random_seq_generator.random_seq_gen(1,P,150)

	#r_top = random_seq_generator.random_seq_gen(1,Prand,139)
	#r_bot= random_seq_generator.random_seq_gen(1,Prand,150)



	## for calculation without GAPDH or fLUC

	#r_top = random_seq_generator.random_seq_gen(1,P,40)
	#r_bot= random_seq_generator.random_seq_gen(1,Prand,98)

	#r_top = random_seq_generator.random_seq_gen(1,Prand,40)
	#r_bot= random_seq_generator.random_seq_gen(1,P,98)

	#r_top = random_seq_generator.random_seq_gen(1,Prand,40)
	#r_bot= random_seq_generator.random_seq_gen(1,Prand,98)


	#r_top = random_seq_generator.random_seq_gen(1,P,58)
	#r_bot= random_seq_generator.random_seq_gen(1,Prand,98)

	#r_top = random_seq_generator.random_seq_gen(1,Prand,58)
	#r_bot= random_seq_generator.random_seq_gen(1,P,98)

	#r_top = random_seq_generator.random_seq_gen(1,Prand,58)
	#r_bot= random_seq_generator.random_seq_gen(1,Prand,98)


	#r_top = random_seq_generator.random_seq_gen(1,P,86)
	#r_bot= random_seq_generator.random_seq_gen(1,Prand,98)

	#r_top = random_seq_generator.random_seq_gen(1,Prand,86)
	#r_bot= random_seq_generator.random_seq_gen(1,P,98)

	r_top = random_seq_generator.random_seq_gen(1,Prand,86)
	r_bot= random_seq_generator.random_seq_gen(1,Prand,98)





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


		## for calculation without GAPDH or fLUC
		#n1 = (float(b[j])/float(98))-(float(t[j])/float(40))
		#n2 = ((r_bot[j][2])**2/float(98)**2)/float(98)
		#n3 = ((r_top[j][2])**2/float(40)**2)/float(40)

		#n1 = (float(b[j])/float(98))-(float(t[j])/float(58))
		#n2 = ((r_bot[j][2])**2/float(98)**2)/float(98)
		#n3 = ((r_top[j][2])**2/float(58)**2)/float(58)

		n1 = (float(b[j])/float(98))-(float(t[j])/float(86))
		n2 = ((r_bot[j][2])**2/float(98)**2)/float(98)
		n3 = ((r_top[j][2])**2/float(86)**2)/float(86)

		
		# to prevent divide by 0
		if (n2+n3) == 0: 
			n4 = 0
		else:
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

