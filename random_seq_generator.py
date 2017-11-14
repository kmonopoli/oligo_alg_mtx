import numpy as np

BASES = ("A","U","C","G")


def random_seq_gen(length,P,size): #P(A),P(U),P(C),P(G)
	total_A = []
	total_U = []
	total_C = []
	total_G = []
	i=0
	while i<size:
		ls = []
		j=0
		while j<size:
			ls.append( ''.join(np.random.choice(BASES, p=P) for _ in range(length)))
			j+=1
		total_A.append(sum([x.count('A') for x in ls]))
		total_U.append(sum([x.count('U') for x in ls]))
		total_C.append(sum([x.count('C') for x in ls]))
		total_G.append(sum([x.count('G') for x in ls]))
		i+=1
	return [[np.mean(total_A),np.median(total_A),np.std(total_A)],[np.mean(total_U),np.median(total_U),np.std(total_U)],[np.mean(total_C),np.median(total_C),np.std(total_C)],[np.mean(total_G),np.median(total_G),np.std(total_G)]]


def make_rand_seqs(length,P,size):
	ls = []
	i=0
	while i<size:
		ls.append( ''.join(np.random.choice(BASES, p=P) for _ in range(length)))
		i+=1
	return ls

MODIFICATIONS = ("f","m","non")

def random_mod_seq_gen(length,P,size): #P(f),P(m),P(non)
	total_f = []
	total_m = []
	total_non = []
	i=0
	while i<size:
		ls = []
		j=0
		while j<size:
			ls.append( ''.join(np.random.choice(MODIFICATIONS, p=P) for _ in range(length)))
			j+=1
		total_f.append(sum([x.count('f') for x in ls]))
		total_m.append(sum([x.count('m') for x in ls]))
		total_non.append(sum([x.count('non') for x in ls]))
		i+=1
	return [[np.mean(total_f),np.median(total_f),np.std(total_f)],[np.mean(total_m),np.median(total_m),np.std(total_m)],[np.mean(total_non),np.median(total_non),np.std(total_non)]]

# Fluoros - always U or C. Me - can be anything. OH - always A or G. 



def random_seq_and_mod_gen(length,P_seq,P_mod_pyr,P_mod_pur,size):
	total_A = []
	total_U = []
	total_C = []
	total_G = []

	total_f_pyr = []
	total_m_pyr = []
	total_non_pyr = []

	total_f_pur = []
	total_m_pur = []
	total_non_pur = []
	
	i=0
	while i<size:
		ls_seq = []
		##ls_mod = []
		ls_mod_pyr = []
		ls_mod_pur = []
		j=0
		while j<size:
			k = 0
			while k<length:
				base = (''.join(np.random.choice(BASES, p=P_seq)))
				if base == 'C' or base == 'U': # CU pyr
					##mod = (''.join(np.random.choice(MODIFICATIONS, p=P_mod_pyr))) 
					mod_pyr = (''.join(np.random.choice(MODIFICATIONS, p=P_mod_pyr))) 
					mod_pur = 'Y'
				else: # AG pur
					##mod = (''.join(np.random.choice(MODIFICATIONS, p=P_mod_pur))) 
					mod_pur = (''.join(np.random.choice(MODIFICATIONS, p=P_mod_pur))) 
					mod_pyr = 'X'
				ls_seq.append(base)
				##ls_mod.append(mod)
				ls_mod_pyr.append(mod_pyr)
				ls_mod_pur.append(mod_pur)
				k+=1
			j+=1
		total_A.append(sum([x.count('A') for x in ls_seq]))
		total_U.append(sum([x.count('U') for x in ls_seq]))
		total_C.append(sum([x.count('C') for x in ls_seq]))
		total_G.append(sum([x.count('G') for x in ls_seq]))

		##total_f.append(sum([x.count('f') for x in ls_mod]))
		##total_m.append(sum([x.count('m') for x in ls_mod]))
		##total_non.append(sum([x.count('non') for x in ls_mod]))


		total_f_pyr.append(sum([x.count('f') for x in ls_mod_pyr]))
		total_m_pyr.append(sum([x.count('m') for x in ls_mod_pyr]))
		total_non_pyr.append(sum([x.count('non') for x in ls_mod_pyr]))


		total_f_pur.append(sum([x.count('f') for x in ls_mod_pur]))
		total_m_pur.append(sum([x.count('m') for x in ls_mod_pur]))
		total_non_pur.append(sum([x.count('non') for x in ls_mod_pur]))
		i+=1



	return [
		[np.mean(total_A),np.median(total_A),np.std(total_A)],
		[np.mean(total_U),np.median(total_U),np.std(total_U)],
		[np.mean(total_C),np.median(total_C),np.std(total_C)],
		[np.mean(total_G),np.median(total_G),np.std(total_G)],

		[np.mean(total_f_pyr),np.median(total_f_pyr),np.std(total_f_pyr)],
		[np.mean(total_m_pyr),np.median(total_m_pyr),np.std(total_m_pyr)],
		[np.mean(total_non_pyr),np.median(total_non_pyr),np.std(total_non_pyr)],

		[np.mean(total_f_pur),np.median(total_f_pur),np.std(total_f_pur)],
		[np.mean(total_m_pur),np.median(total_m_pur),np.std(total_m_pur)],
		[np.mean(total_non_pur),np.median(total_non_pur),np.std(total_non_pur)]]


PURINE_PYRIM_MODS = ("f_pyr","m_pyr","non_pyr","f_pur","m_pur","non_pur")
def random_mod_gen_pyrim_purine(length,P,size):
	total_f_pyr = []
	total_m_pyr = []
	total_non_pyr = []

	total_f_pur = []
	total_m_pur = []
	total_non_pur = []
	
	i=0
	while i<size:
		ls_mod = []
		j=0
		while j<size:
			k = 0
			while k<length:
				mod_pyrim_pur = ''.join(np.random.choice(PURINE_PYRIM_MODS, p=P))
				ls_mod.append(mod_pyrim_pur)
				k+=1
			j+=1

		total_f_pyr.append(sum([x.count('f_pyr') for x in ls_mod]))
		total_m_pyr.append(sum([x.count('m_pyr') for x in ls_mod]))
		total_non_pyr.append(sum([x.count('non_pyr') for x in ls_mod]))


		total_f_pur.append(sum([x.count('f_pur') for x in ls_mod]))
		total_m_pur.append(sum([x.count('m_pur') for x in ls_mod]))
		total_non_pur.append(sum([x.count('non_pur') for x in ls_mod]))
		i+=1



	return [
		[np.mean(total_f_pyr),np.median(total_f_pyr),np.std(total_f_pyr)],
		[np.mean(total_m_pyr),np.median(total_m_pyr),np.std(total_m_pyr)],
		[np.mean(total_non_pyr),np.median(total_non_pyr),np.std(total_non_pyr)],

		[np.mean(total_f_pur),np.median(total_f_pur),np.std(total_f_pur)],
		[np.mean(total_m_pur),np.median(total_m_pur),np.std(total_m_pur)],
		[np.mean(total_non_pur),np.median(total_non_pur),np.std(total_non_pur)]
		]


PER_BASE_MODS = ("f_A","m_A","non_A","f_U","m_U","non_U","f_C","m_C","non_C","f_G","m_G","non_G")
def random_mod_gen_per_base(length,P,size):
	total_f_A = []
	total_m_A = []
	total_non_A = []

	total_f_U = []
	total_m_U = []
	total_non_U = []
	
	total_f_C = []
	total_m_C = []
	total_non_C = []
	
	total_f_G = []
	total_m_G = []
	total_non_G = []

	
	i=0
	while i<size:
		ls_mod = []
		j=0
		while j<size:
			k = 0
			while k<length:
				mod_perBase = ''.join(np.random.choice(PER_BASE_MODS, p=P))
				ls_mod.append(mod_perBase)
				k+=1
			j+=1

		total_f_A.append(sum([x.count('f_A') for x in ls_mod]))
		total_m_A.append(sum([x.count('m_A') for x in ls_mod]))
		total_non_A.append(sum([x.count('non_A') for x in ls_mod]))

		total_f_U.append(sum([x.count('f_U') for x in ls_mod]))
		total_m_U.append(sum([x.count('m_U') for x in ls_mod]))
		total_non_U.append(sum([x.count('non_U') for x in ls_mod]))
		
		total_f_C.append(sum([x.count('f_C') for x in ls_mod]))
		total_m_C.append(sum([x.count('m_C') for x in ls_mod]))
		total_non_C.append(sum([x.count('non_C') for x in ls_mod]))
		
		total_f_G.append(sum([x.count('f_G') for x in ls_mod]))
		total_m_G.append(sum([x.count('m_G') for x in ls_mod]))
		total_non_G.append(sum([x.count('non_G') for x in ls_mod]))
		i+=1



	return [
		[np.mean(total_f_A),np.median(total_f_A),np.std(total_f_A)],
		[np.mean(total_m_A),np.median(total_m_A),np.std(total_m_A)],
		[np.mean(total_non_A),np.median(total_non_A),np.std(total_non_A)],

		[np.mean(total_f_U),np.median(total_f_U),np.std(total_f_U)],
		[np.mean(total_m_U),np.median(total_m_U),np.std(total_m_U)],
		[np.mean(total_non_U),np.median(total_non_U),np.std(total_non_U)],

		[np.mean(total_f_C),np.median(total_f_C),np.std(total_f_C)],
		[np.mean(total_m_C),np.median(total_m_C),np.std(total_m_C)],
		[np.mean(total_non_C),np.median(total_non_C),np.std(total_non_C)],

		[np.mean(total_f_G),np.median(total_f_G),np.std(total_f_G)],
		[np.mean(total_m_G),np.median(total_m_G),np.std(total_m_G)],
		[np.mean(total_non_G),np.median(total_non_G),np.std(total_non_G)]
		]

## Testing ##

#print random_seq_and_mod_gen(1,[0.5,0.0,0.0,0.5],[1.0,0.0,0.0],[0.25,0.25,0.5],10)










