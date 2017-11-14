## Data type for Duplex data


class Duplex:
	iD = 0
	gene = ''
	datePrep = 'mm.dd.yyyy'
	twntyMer = ''
	eightyMer = ''
	fm = bool
	expP = 0.0
	stdV = 0.0
	posCntP = 0.0
	ntcP = 0.0
	dateExper = 'mm.dd.yyyy'
	cellLine = ''
	screenType = ''
	sense_fm_seq = ''
	antisense_fm_seq = ''
	gene_region = '' # 45mer gene region
	twntyMerPosnMtrx = []
	
	def __init__(self, iD, duplex_name, twntyMer, fm_val, eightyMer, datePrep, expP, stdV, posCntP, ntcP, dateExper, cellLine, screenType, sense_fm_seq, antisense_fm_seq):
		self.iD = int(iD)
		self.gene = self.get_geneName(duplex_name)
		self.datePrep = self.get_date(datePrep)
		self.twntyMer = twntyMer
		self.eightyMer = eightyMer
		self.fm = self.get_fm(fm_val)
		self.expP = float(expP)
		self.stdV = float(stdV)
		self.ntcP = float(ntcP)
		self.dateExper = self.get_date(dateExper)
		self.cellLine = cellLine
		self.screenType = screenType
		self.posCntP = self.get_posCntP(posCntP)
		self.sense_fm_seq = self.make_mod_ls(sense_fm_seq)	
		self.antisense_fm_seq = self.make_mod_ls(antisense_fm_seq)
		self.set_gene_region()

	def get_geneName(self,duplex_name):
		return duplex_name[:duplex_name.find('_')]
				
	def get_fm(self,fm_val):
		if fm_val == 'fm':
			return True
		elif fm_val == '':
			return False
		else:
	
			print "ERROR: for Duplex ID: ", self.iD
			print "fm value from file is ",fm_val," not a blank or fm"
			quit()
	
	def get_date(self,date):
		if len(date) != len('mm.dd.yyyy'):
			print "ERROR: for Duplex ID: ", self.iD
			print "Problem in date ",date," not correct formatting based on length of string"
			quit()
		else:
			return date
		#try:
		#	int(date.replace('.',','))
		#except:
		#	print "ERROR: date ",date," not correct formatting when removed '.' must be an int"
		#	quit()
	
	def get_posCntP(self,posCntP):
		if self.screenType == 'qPCR':
			self.posCntP = None
		elif self.screenType == 'Renilla-Firefly Luciferase':
			try: 
				self.posCntP = float(posCntP)
			except:
		#		print "WARNING: for Duplex ID: ", self.iD
		#		print "Postive Control % set to None for Luciferase screen"
				self.posCntP = None
		else:
			print "ERROR: for Duplex ID: ", self.iD
			print "Could not get Positive Control % because Screen Type ", self.screenType , " is not of type qPCR or Renilla-Firefly Luciferase"
			quit()

	def get_data(self):
		print 'ID: ', self.iD
		print 'gene: ',self.gene
		print 'Date prepared: ', self.datePrep
		print '20mer: ', self.twntyMer
		print '80mer: ', self.eightyMer
		print 'fm: ',self.fm
		print '% Expr: ', self.expP
		print 'Stdev: ', self.stdV
		print 'PosCnt%: ', self.posCntP
		print 'NTC%: ', self.ntcP
		print 'Date experiment: ',self.dateExper
		print 'Cell line: ', self.cellLine
		print 'Screen Type: ',self.screenType
		print 'fm seq: ',self.fm_seq
		print 'gene region: ',self.gene_region
		return
	
	def set_gene_region(self):
		self.gene_region = self.eightyMer[15:60]
		return
	
	def make_mod_ls(self,seq):
		seq = seq.replace('\n','')
		seq = seq.replace('#','.')
		return seq.split('.')

## Testing
'''
import csv_io

x = csv_io.CSV_Reader('output_data_diana_less_than_80_removed.csv')#'output_data.csv')
n= x.get_Duplex_data(100)

d = Duplex(*n)


print d.gene
print d.iD
print d.expP
print d.twntyMer
print d.gene_region
print len(d.gene_region)
'''


