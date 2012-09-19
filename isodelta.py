#!/usr/bin/python

import os

##############################################################################
##### Extract results from a given folder ####################################
##############################################################################
def extractParameter (filename, name):
	parafile=open(filename, 'r')
	lines=parafile.readlines()
	parafile.close()
	for l in lines:
		if l.startswith(name):
			value=l.split()[2]
			return value
	print "Warning: parameter %s in file %s not found!" % (name, filename)
	exit(1)

def extractResultValue (filename):
	f=open(filename, 'r')
	l=f.readline()
	f.close()
	return l

def isResults (resultsFolder):
	parafilename="%s/parameter.cfg" % resultsFolder
	mufilename="%s/results/mu.dat" % resultsFolder
	if os.path.exists(parafilename) and os.path.exists(mufilename):
		return True
	else:
		return False
def extractResults (resultsFolder):
	parafilename="%s/parameter.cfg" % resultsFolder
	mufilename="%s/results/mu.dat" % resultsFolder
	mu=float(extractResultValue(mufilename))
	Delta=-mu
	nc=float(extractParameter(parafilename, 'concentration'))
	T=float(extractParameter(parafilename, 'temperature'))
	N=int(extractParameter(parafilename, 'N'))
	impurity=extractParameter(parafilename, 'impurity')
	material=''
	if impurity=='Gd':
		material='EuO'
	elif impurity=='None':
		material='Sub'
	else:
		print "Error: Unable to extract material type from folder %s" % resultsFolder
		exit(1)
		
	parafilename="%s/parameter.cfg" % resultsFolder
	return (material, N, nc, T, Delta);

##############################################################################
##### Database for energy shifts in EuO and substrate ########################
##############################################################################
def string(val):
	if type(val).__name__=='str':
		return val
	elif type(val).__name__=='int':
		return "%d" % val
	elif type(val).__name__=='float':
		return "%0.15e" % val
	else:
		print "Error: Sting conversion: Unknown type. Break"
		exit(1)
class database:
	def __init__(self):
		self.names=('material', 'N', 'nc', 'T', 'Delta')

	def write(self, filename):
		f=open(filename, 'w')
		f.write('#mat\tN\tnc\t\t\tT\t\t\tDelta\n')

		for d in self.data:
			for val in d:
				f.write(string(val))
				f.write('\t')
			f.write('\n')

	def set(self, d):
		self.data=d
		self.names=('material', 'N', 'nc', 'T', 'Delta')

	def read(self, filename):
		f=open(filename, 'r')
		lines=f.readlines()
		f.close()

		self.data=[]
		for l in lines[2:]:
			d=l.split()
			material_val=d[0]
			N_val=int(d[1])
			nc_val=float(d[2])
			T_val=float(d[3])
			Delta_val=float(d[4])
			self.data.append((material_val, N_val, nc_val, T_val, Delta_val))
	
	def exists(self, material, N, nc, T):
		for d in self.data:
			if material==d[0] and N==d[1] and nc==d[2] and T==d[3]:
				return True
		return False

	def getDelta(self, material, N, nc, T):
		for d in self.data:
			if material==d[0] and N==d[1] and nc==d[2] and T==d[3]:
				return d[4]

		print "Error: Unable to find dataset in database:"
		print "material=%s" % material
		print "N=%s" % N
		print "nc=%s" % nc 
		print "T=%s" % T
		print "Break."
		exit(1)
	
	def extract(self, topResultsFolders):
		self.data=[]
		for topfolder in topResultsFolders:
			for d in os.listdir(topfolder):
				folder=os.path.join(topfolder, d)
				if os.path.isdir(folder) and isResults(folder):
					self.data.append(extractResults(folder))

		# sort by temperature
		self.data=sorted(self.data, key = lambda element : element[0])
