## Datatypes for Reading and Writing to CSVs
import csv

class CSV_Reader:
	
	file_name = ''
	file_data_ls_unparsed = []
	headers = []
	file_data = [] # no headers
	lines = 0 # not including header
	columns = 0
	
	def __init__(self,file_name):
		if file_name[-4:] == '.csv':
			self.file_name = file_name
		else:
			self.file_name = file_name+'.csv'

		self.file_data_ls_unparsed = self.read_file()
		self.headers = [x for x in self.file_data_ls_unparsed[0] if x != '']
		self.file_data = self.file_data_ls_unparsed[1:]
		self.lines = len(self.file_data)				
		self.columns = len(self.headers)

	def read_file(self):
		data = []
		with open(self.file_name,'rb') as f:
			data = f.readlines()
			data = [x.replace('\r\n',',').split(',') for x in data]
			data = [x for x in data]
			f.close()
		return data
		
	def get_row(self,n):
		if n > self.lines:
			print "ERROR: get_row request for line ",n," with only ",self.lines," lines in file"
			quit() 
		else:
			return self.file_data[n]

	def get_columns(self):
		ls = [[]]*self.columns
		ls = [[x[i] for x in  self.file_data] for i in range(self.columns)]
		return ls
	
	def get_Duplex_data(self,n):
		data = []
		row = self.get_row(n)
		index = [0,2,8,10,11,13,14,15,16,17,18,19,20,22,23] # list of indexes in row with desired Duplex data
		for i in index:
			data.append(row[i])
		return data

class CSV_Writer:
	
	file_name = ''
	headers = []
	file_data = [] # no headers
	lines = 0 # not including header
	columns = 0
	rows = 0
	inRows = bool # true if data_ls passed to constructor in rows, false if data passed is in collumns
	
	def __init__(self,file_name,headers,data_ls,inRows): # data_ls is list of rows or columns 
		if file_name[-4:] == '.csv':
			self.file_name = file_name
		else:
			self.file_name = file_name+'.csv'

		self.headers = headers
		self.file_data = data_ls
		self.inRows = inRows

	def write_to_file(self):
		if self.inRows:
			self.write_rows()
		else:
			self.write_cols()	

	def write_rows(self):
		with open(self.file_name, 'wb') as csvfile:
			fwriter = csv.writer(csvfile, delimiter=',')
			fwriter.writerow(self.headers)
			for x in self.file_data:
				fwriter.writerow(x)
			print "Data written to ",self.file_name
		return
	
	def write_cols(self):
		with open(self.file_name, 'wb') as csvfile:
			fwriter = csv.writer(csvfile, delimiter=',')
			fwriter.writerow(self.headers)
			ls = zip(*self.file_data)
			ls = [list(x) for x in ls]
			for x in ls:
				fwriter.writerow(x)
			print "Data written to ",self.file_name
		return

		


## Testing
'''
x = CSV_Reader('qPCR_nonfm_normalized_data.csv')

ls1 =  x.get_row(3)#x.get_Duplex_data(3)
ls2 = x.headers

print x.get_Duplex_data(3)
i = 0
while i<len(ls2):
	print i," : ",ls1[i]," : ",ls2[i]
	i+=1
'''


