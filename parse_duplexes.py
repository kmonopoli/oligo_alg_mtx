#!/usr/bin/python
## Analyzes and sorts Duplex data

import csv_io
import duplex

class Duplex_Analysis:
	all_duplexes = []

	luc = []
	qpcr = []

	fm = []
	Nfm = []

	luc_fm = []
	luc_Nfm = []
	qpcr_fm = []
	qpcr_Nfm = []
	
	# can only hold one assortment at a time!
	hit_bin = []
	fail_bin = []

	def __init__(self, fName):
		x = csv_io.CSV_Reader(fName)
		ls = [ x.get_Duplex_data(n) for n in range(x.lines) ]
		self.all_duplexes = [ duplex.Duplex(*n) for n in ls ]

		self.luc = self.get_by_screen( self.all_duplexes, 'Renilla-Firefly Luciferase' )
		self.qpcr = self.get_by_screen( self.all_duplexes, 'qPCR' )
		self.fm = self.get_by_fm( self.all_duplexes, True )
		self.Nfm = self.get_by_fm( self.all_duplexes, False )				

		self.luc_fm = self.get_by_screen( self.fm, 'Renilla-Firefly Luciferase')
		self.luc_Nfm = self.get_by_screen( self.Nfm, 'Renilla-Firefly Luciferase')
		self.qpcr_fm = self.get_by_screen( self.fm, 'qPCR')
		self.qpcr_Nfm = self.get_by_screen( self.Nfm, 'qPCR')


	def get_by_screen(self, ls, screen):
		data = [n for n in ls if n.screenType == screen]
		return data
	
	def get_by_fm(self, ls, fm_bool):
		data = [n for n in ls if n.fm == fm_bool]
		return data
	
	def sort_by_expP(self, ls):
		# returns list from lowest Expression % to highest
		sorted_ls = sorted(ls, key=lambda x: x.expP, reverse=False)		
		return sorted_ls

	def bin_by_expP(self, ls, p_threshold_hit, p_threshold_fail):
		# a method Diana used, separates into bins of values with expressions above and below given values
		hits = []
		fails = []
		ls = self.sort_by_expP(ls)
		hits = [x for x in ls if x.expP <= p_threshold_hit]
		fails = [x for x in ls if x.expP >= p_threshold_fail]
		self.hit_bin = hits
		self.fail_bin = fails
		return [hits,fails] 

	def bin_by_number(self, ls, num_hit, num_fail): # hit is low Expression %
		# other method Diana used; bins into two given sized groups of hits and fails
		if len(ls) < num_hit+num_fail:
			print "ERROR: cannot bin, not enough data"
			quit()
		hits = []
		fails = []
		ls = self.sort_by_expP(ls)
		hits = ls[:num_hit]
		fails = ls[-num_fail:]
		self.hit_bin = hits
		self.fail_bin = fails
		return [hits,fails] 

	def get_seqNoMods(self,ls):
		newLs = []
		for x in ls:
			n = ""
			for y in x:

				q = y.replace("f","").replace("m","").replace("P","").replace("#","")
				n=n+q
			newLs.append(n)	
		return newLs


	def bin_by_number_unsorted(self, ls, num_hit, num_fail): # hit is low Expression %
		# for fluch random pvalue generator
		if len(ls) < num_hit+num_fail:
			print "ERROR: cannot bin, not enough data"
			quit()
		hits = []
		fails = []
		#ls = self.sort_by_expP(ls)
		hits = ls[:num_hit]
		fails = ls[-num_fail:]
		self.hit_bin = hits
		self.fail_bin = fails
		return [hits,fails] 

	def bin_by_rankOrder(self,ls,hit_threshold,fail_threshold):
		# takes all data and bins by rank order FINISH
		ls = self.bin_by_experiment(ls)
		new_ls = []
		i = 0
		while i<len(ls):
			new_ls.append(sorted(ls[i], key=lambda x: x.expP, reverse=False))
			i+=1
		return new_ls

	def bin_by_experiment(self,ls):
		new_ls = [] # holds binned Duplexes
		dates = [] # keeps track of positions of dates added in new_ls
		genes = [] # keeps track of positions of genes added in new_ls
		i = 0
		while i < len(ls):
			d = ls[i]		
			indices_1 = [n for n, x in enumerate(dates) if x == d.dateExper]
			indices_2 = [n for n, x in enumerate(genes) if x == d.gene]
			intersect = [x for x in indices_1 if x in indices_2]
			if len(intersect)==1:
				new_ls[intersect[0]].append(d)
			elif len(intersect)==0:
				new_ls.append([d])
				dates.append(d.dateExper)
				genes.append(d.gene)
			else:
				print "ERROR: more than one experiment with date: ",d.dateExper," and gene: ",d.gene
				quit()
			i+=1
		if len(ls) != sum([len(n) for n in new_ls]):
			print "ERROR: data was lost, staring amount was ",len(ls)," Duplexes, ended with ",sum([len(n) for n in ls])
			quit()
		return new_ls
		

	def get_freqsPerPos(self, ls):
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
	
	def normalize_seq_data(self,ls):
		s= ls[0][0]+ls[1][0]+ls[2][0]+ls[3][0]
		return [[float(n)/s for n in k] for k in ls]	

	def subtract_seq_data(self,hits,fails):
		# normalize first   
		subtracted = []
		k = 0
		while k < len(hits):
			l = [j-i for i, j in zip(hits[k],fails[k])]
			subtracted.append(l)
			k+=1
		return subtracted


	def get_mod_freqsPerPos(self, ls):
		freqf = [0]*len(ls[0])
		freqm = [0]*len(ls[0])
		freqnon = [0]*len(ls[0])
		for seq in ls:
			i = 0
			while i < len(seq):
				s = seq[i].replace('G','').replace('A','').replace('U','').replace('C','').replace('1','').replace('P','')

				if s == 'f':
					freqf[i]+=1	
		
				elif s == 'm':
					freqm[i]+=1	
		
				elif s == '':
					freqnon[i]+=1
				else:
					print "ERROR: sequence ", seq ," contains non f/m modification: ",s	
				i+=1
		return [freqf,freqm,freqnon]

				
	def normalize_mod_seq_data(self,ls):
		s= ls[0][0]+ls[1][0]+ls[2][0]
		return [[float(n)/s for n in k] for k in ls]	
		

	def make_20mer_antisense(self,ls):
		## takes reverse complement of 20mer targetting region and returns sequence of antisense oligo
		ls = [x.replace("G","c").replace("C","g").replace("A","u").replace("U","a").upper()	for x in ls]
		ls = [x1[::-1] for x1 in ls]
		## because every antisense oligo starts with U:
		ls = [x2[1:] for x2 in ls]
		ls = ['U'+x3 for x3 in ls]
		return ls



	def get_mod_freqsPerPos_purines_vs_pyrimidines(self, ls):
		# pyrimidines  CU
		freqf_pyr = [0]*len(ls[0])
		freqm_pyr = [0]*len(ls[0])
		freqnon_pyr = [0]*len(ls[0])
		# purines AG
		freqf_pur = [0]*len(ls[0]) # will always be zeros
		freqm_pur = [0]*len(ls[0])
		freqnon_pur = [0]*len(ls[0])

		for seq in ls:
			i = 0
			while i < len(seq):
				s = seq[i].replace('1','').replace('P','') ## may have to account for these mods later!
				if 'C' in s or 'U' in s:
					s = s.replace('C','').replace('U','')
					if s == 'f':
						freqf_pyr[i]+=1	
					elif s == 'm':
						freqm_pyr[i]+=1	
					elif s == '':
						freqnon_pyr[i]+=1
					else:
						print "ERROR: sequence ", seq ," contains non f/m modification: ",s	
				elif 'G' in s or 'A' in s:
					s = s.replace('G','').replace('A','')
					if s == 'm':
						freqm_pur[i]+=1	
					elif s == '':
						freqnon_pur[i]+=1
					else:
						print "ERROR: sequence ", seq ," contains non f/m modification: ",s	
				else:
					print "ERROR: sequence ",seq ," contains non AUGC bases: ", seq[i]
				i+=1
		return [[freqf_pyr,freqm_pyr,freqnon_pyr],[freqf_pur,freqm_pur,freqnon_pur]]

	def normalize_mod_seq_data_purine_pyrimidine(self,ls):
		p_ls = [[[],[],[]],[[],[],[]]]
		i = 0
		while i < len(ls): # 2
			j = 0
			while j < len(ls[i][0]): # 20
				s = ls[i][0][j]+ls[i][1][j]+ls[i][2][j]
				if s == 0:
					p_ls[i][0].append(0.0);	
					p_ls[i][1].append(0.0);	
					p_ls[i][2].append(0.0);
				else:
					p_ls[i][0].append(float(ls[i][0][j])/float(s))	
					p_ls[i][1].append(float(ls[i][1][j])/float(s))	
					p_ls[i][2].append(float(ls[i][2][j])/float(s))	
				j+=1
			i+=1	
		return p_ls		


	def get_mod_freqsPerPos_purines_vs_pyrimidines_2(self, ls):
		# pyrimidines  CU
		freqf_pyr = [0]*len(ls[0])
		freqm_pyr = [0]*len(ls[0])
		freqnon_pyr = [0]*len(ls[0])
		# purines AG
		freqf_pur = [0]*len(ls[0]) # will always be zeros
		freqm_pur = [0]*len(ls[0])
		freqnon_pur = [0]*len(ls[0])

		for seq in ls:
			i = 0
			while i < len(seq):
				s = seq[i].replace('1','').replace('P','') ## may have to account for these mods later!
				if 'C' in s or 'U' in s:
					s = s.replace('C','').replace('U','')
					if s == 'f':
						freqf_pyr[i]+=1	
					elif s == 'm':
						freqm_pyr[i]+=1	
					elif s == '':
						freqnon_pyr[i]+=1
					else:
						print "ERROR: sequence ", seq ," contains non f/m modification: ",s	
				elif 'G' in s or 'A' in s:
					s = s.replace('G','').replace('A','')
					if s == 'm':
						freqm_pur[i]+=1	
					elif s == '':
						freqnon_pur[i]+=1
					else:
						print "ERROR: sequence ", seq ," contains non f/m modification: ",s	
				else:
					print "ERROR: sequence ",seq ," contains non AUGC bases: ", seq[i]
				i+=1
		return [freqf_pyr,freqm_pyr,freqnon_pyr,freqf_pur,freqm_pur,freqnon_pur]


	def normalize_mod_seq_data_purine_pyrimidine_2(self,ls):
		s= ls[0][0]+ls[1][0]+ls[2][0]+ls[3][0]+ls[4][0]+ls[5][0]
		return [[float(n)/s for n in k] for k in ls]	


	def get_mod_freqsPerPos_per_base(self, ls):
		# pyrimidines  A
		freqf_A = [0]*len(ls[0])
		freqm_A = [0]*len(ls[0])
		freqnon_A = [0]*len(ls[0])
		# pyrimidines  U
		freqf_U = [0]*len(ls[0])
		freqm_U = [0]*len(ls[0])
		freqnon_U = [0]*len(ls[0])
		# pyrimidines  C
		freqf_C = [0]*len(ls[0])
		freqm_C = [0]*len(ls[0])
		freqnon_C = [0]*len(ls[0])
		# pyrimidines  G
		freqf_G = [0]*len(ls[0])
		freqm_G = [0]*len(ls[0])
		freqnon_G = [0]*len(ls[0])



		for seq in ls:
			i = 0
			while i < len(seq):
				s = seq[i].replace('1','').replace('P','') ## may have to account for these mods later!
				if 'A' in s:
					s = s.replace('A','')
					if s == 'f':
						freqf_A[i]+=1	
					elif s == 'm':
						freqm_A[i]+=1	
					elif s == '':
						freqnon_A[i]+=1
					else:
						print "ERROR: sequence ", seq ," contains non f/m modification: ",s	
				elif 'U' in s:
					s = s.replace('U','')
					if s == 'f':
						freqf_U[i]+=1	
					elif s == 'm':
						freqm_U[i]+=1	
					elif s == '':
						freqnon_U[i]+=1
					else:
						print "ERROR: sequence ", seq ," contains non f/m modification: ",s	
				elif 'C' in s:
					s = s.replace('C','')
					if s == 'f':
						freqf_C[i]+=1	
					elif s == 'm':
						freqm_C[i]+=1	
					elif s == '':
						freqnon_C[i]+=1
					else:
						print "ERROR: sequence ", seq ," contains non f/m modification: ",s	
				elif 'G' in s:
					s = s.replace('G','')
					if s == 'f':
						freqf_G[i]+=1	
					elif s == 'm':
						freqm_G[i]+=1	
					elif s == '':
						freqnon_G[i]+=1
					else:
						print "ERROR: sequence ", seq ," contains non f/m modification: ",s	

				else:
					print "ERROR: sequence ",seq ," contains non AUGC bases: ", seq[i]
				i+=1
		return [
				freqf_A,freqm_A,freqnon_A,
				freqf_U,freqm_U,freqnon_U,
				freqf_C,freqm_C,freqnon_C,
				freqf_G,freqm_G,freqnon_G]

	def normalize_mod_seq_data_per_base(self,ls):
		s= ls[0][0]+ls[1][0]+ls[2][0]+ls[3][0]+ls[4][0]+ls[5][0]+ls[6][0]+ls[7][0]+ls[8][0]+ls[9][0]+ls[10][0]+ls[11][0]
		return [[float(n)/s for n in k] for k in ls]	

### testing
