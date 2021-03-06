#!/usr/bin/python

import os
import subprocess

##############################################################################
##############################################################################
##### Extract results from a given folder ####################################
##############################################################################
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
		# calculation of N layer in mirror symmetric mode
		# corresponds to 2*N-1 layer in non mirror symmetric
		# mode
		N=2*N-1
	else:
		print "Error: Unable to extract material type from folder %s" % resultsFolder
		exit(1)
		
	parafilename="%s/parameter.cfg" % resultsFolder
	return (material, N, nc, T, Delta);

##############################################################################
##############################################################################
##### Database for energy shifts in isolated subsystems EuO and substrate ####
##############################################################################
##############################################################################
# convert data values to strings
def tostring(val):
	if type(val).__name__=='str':
		return val
	elif type(val).__name__=='int':
		return "%d" % val
	elif type(val).__name__=='float':
		return "%0.15e" % val
	else:
		print "Error: Sting conversion: Unknown type. Break"
		exit(1)

class isodeltabase:
	def __init__(self):
		self.names=('material', 'N', 'nc', 'T', 'Delta')

	# write database to file
	def write(self, filename):
		f=open(filename, 'w')
		f.write('#mat\tN\tnc\t\t\tT\t\t\tDelta\n')

		for d in self.data:
			for val in d:
				f.write(tostring(val))
				f.write('\t')
			f.write('\n')

	# fill database 
	def set(self, d):
		self.data=d
		self.names=('material', 'N', 'nc', 'T', 'Delta')

	# read database from file
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

	# read in database from remote file
	def download(self, remotepath='stollenw@stgeorgenamreith.th.physik.uni-bonn.de:/home/stollenw/projects/euo/database/isodelta.db'):
		cmd='scp %s isodelta.db' % remotepath
		try:
			subprocess.call(cmd, shell=True)
			self.read('isodelta.db')
		except:
			print 'Error: Failed to retrieve remote isodelta database: %s' % remotepath
			print 'Break.'
			exit(1)
		
		
	# check if special dataset exists in database
	def exists(self, material, N, nc, T):
		for d in self.data:
			if material==d[0] and N==d[1] and nc==d[2] and T==d[3]:
				return True
		return False

	# extract Delta
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
	
	# fill database by extracting results form a special folder
	def fill(self, topResultsFolders):
		self.data=[]
		for topfolder in topResultsFolders:
			for d in os.listdir(topfolder):
				folder=os.path.join(topfolder, d)
				if os.path.isdir(folder) and isResults(folder):
					self.data.append(extractResults(folder))

		# sort by temperature
		self.data=sorted(self.data, key = lambda element : element[0])

##############################################################################
##############################################################################
##### Class managing the EuO program run commands ############################
##############################################################################
##############################################################################
class run:
	def __init__(self, c):
		self.cmd=c

	######################################################################
	##### Interface run - Add energy shifts of isolated materials ########
	######################################################################
	def get_option(self, key, default=None):
		# get command line options
		options=self.cmd.split()
		# extract temperature (if none given, default is 40 K)
		try:
			i=options.index(key)
			value=options[i+1].rstrip('/')
			return value
		except ValueError:
			if default!=None:
				return default
			else:
				print 'Error: Failed to extract value for option %s in run command: %s. No default was given. Break.' % (key, self.cmd)
				exit(1)
	def option_exists(self, key):
		# get command line options
		options=self.cmd.split()
		# extract temperature (if none given, default is 40 K)
		try:
			i=options.index(key)
			value=options[i+1].rstrip('/')
			return True
		except ValueError:
			return False

	def add_isodeltas(self, dbpath=None):
		# get relevant command line options
		T=float(self.get_option('-t'))
		N=int(self.get_option('-n'))
		N0=int(self.get_option('-m'))
		ni=float(self.get_option('-x', 0.01))
		deltaW=float(self.get_option('--Delta_W', 0.0))
		ncr=float(self.get_option('--n_cr', 0.01))
		Delta_g=float(self.get_option('--Delta_g', 0.125))
		insulator=self.option_exists('--insulator')

		# read database
		db=isodeltabase()
		if dbpath==None:
			db.download()
		else:
			db.read(dbpath)

		# get Deltas (isolated energy shifts) for EuO and metallic substrate
		# for the insulating case, just shift Delta_r by Delta_g
		Delta_l=db.getDelta('EuO', N, ni, T)
		self.cmd+=" --Delta_l0 %0.15e" % Delta_l
		if not insulator:
			Delta_r=db.getDelta('Sub', N0, ncr, T)
			self.cmd+=" --Delta_r0 %0.15e" % Delta_r

	#######################################################################
	##### Add suitable input to euo run command if available ##############
	#######################################################################
	# Add input path to EuO run command 'runcmd'. Hereby the output path and the 
	# temperature in 'runcmd' is extracted. It's parent directory is then searched 
	# for folder containing results with smaller temperatures. The full run command
	# with input options is then returned.
	def add_input (self, path=None):
		options=self.cmd.split()
		# extract output folder from run command if path variable is not given
		if path==None:
			try:
				i=options.index('-o')
				path=os.path.dirname(os.path.abspath(options[i+1].rstrip('/')))
			except ValueError:
				path=os.getcwd()
	
		# extract temperature (if none given, default is 40 K)
		T=0.0
		try:
			i=options.index('-t')
			T=float(options[i+1].rstrip('/'))
		except ValueError:
			T=40.0
	
		inputOptions=''
		# Search search 'resultsFolder' for sub folders containing results with
		# smaller temperatures thant 'T' 
		resultFolders=[]
		for d in os.listdir(path):
			folder=os.path.join(path, d)
			if os.path.isdir(folder) and isResults(folder):
				resultFolders.append((folder, extractResults(folder)[3]))
	
		# find a folder with lower temperature than T
		tmax=0.0
		for (f,t) in resultFolders:
			if t>tmax and t<=T:
				tmax=t
				inputFolder=f
	
		if tmax>0.0:
			self.cmd+=" -i " + inputFolder

db=isodeltabase()
db.fill(['../../runs/runs_version-a177549/pure_ncr0.50/n3/','../../runs/runs_version-a177549/n5/'])
db.write('../../database/isodelta.db')
r=run('euo.out -n 5 -m 5 --n_cr 0.5 -o ../../runs/runs_version-a177549/n5/output_t050/ -t 60')
r.add_input('../../runs/runs_version-a177549/n5_m2_ncr0.10/')
r.add_isodeltas()
print r.cmd
