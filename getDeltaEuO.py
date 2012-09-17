#!/usr/bin/python

import subprocess
from sys import *
import os
from math import *

# Get calculated mu for bulk EuO and return resulting effective
# Band center Delta such, that chemical potential is zero
# Delta = former Delta (=1.0) -mu
def getDeltaEuO (resultsFolder):
	# Delta0 is always 1.0
	Delta0=1.0
	# get mu 
	mufilename="%s/results/mu.dat" % resultsFolder
	if os.path.exists(mufilename):
		mufile=open(mufilename, 'r')
		l=mufile.readline()
		mufile.close()
		mu=float(l)
		Delta=Delta0-mu
	else:
		print "Error: getDeltaEuO: Path '%s' does not exist. Break." % mufilename 
		#exit(1)
		return ( -1,-1,-1);
	# get temperature 
	T=-1.0
	ni=-1.0
	foundT=False
	foundNi=False
	parafilename="%s/parameter.cfg" % resultsFolder
	if os.path.exists(parafilename):
		parafile=open(parafilename, 'r')
		lines=parafile.readlines()
		parafile.close()
		for l in lines:
			if l.startswith('temperature'):
				T=float(l.split()[2])
				foundT=True
			if l.startswith('concentration'):
				ni=float(l.split()[2])
				foundNi=True
		if not foundT:
			print "Error: getDeltaEuO: Could not extract temperature. Break." % mufilename 
			exit(1)
		if not foundNi:
			print "Error: getDeltaEuO: Could not extract impurity concentration. Break." % mufilename 
			exit(1)
	else:
		print "Error: getDeltaEuO: Path '%s' does not exist. Break." % mufilename 
		exit(1)
	return ( 1.0 - mu, T, ni);
	
#main program
def main():
	if (len(argv))<1+2:
		stderr.write( """not enought arguments. call with
		1.) Output filename 
		2.) Result folder of the bulk runs
		\n""")
		exit(1)
	#parameter
	ofn=argv[1]
	topResultsFolder=argv[2]

	euodata=[]
	for d in os.listdir(topResultsFolder):
		if os.path.isdir(os.path.join(topResultsFolder, d)):
			euodata.append( getDeltaEuO(os.path.join(topResultsFolder, d)))

	# sort by temperature
	euodata=sorted(euodata, key = lambda element : element[1])
	# add Delta for substrate and write into file
	of=open(ofn, 'w')
	for data in euodata:
		(Delta_l, T, ni)=data
		process=subprocess.Popen('./freedos.out -t %f -n %f' % (T, ni), shell=True, stdout=subprocess.PIPE)
		Delta_r=float(process.communicate()[0].rstrip('\n'))
		print T, ni, Delta_l, Delta_r
		of.write('%0.15e\t%0.15e\t%0.15e\t%0.15e\n' % (T, ni, Delta_l, Delta_r))
	of.close()

if __name__=="__main__":
	main()
