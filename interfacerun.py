#!/usr/bin/python

import isodelta
import argparse

# Update database for energy shifts in EuO and substrate
def main():
	parser = argparse.ArgumentParser(description='Update database for energy shifts in EuO and substrate')
  	parser.add_argument('-d', '--database', type=str, help='Database file name', default='isodelta.dat')
  	parser.add_argument('-n', '--N', type=int, help='Number of EuO layers', default=1)
  	parser.add_argument('-m', '--N0', type=int, help='Number of Substrate layers', default=0)
  	parser.add_argument('-o', '--output', type=str, help='Output folder')

	args = parser.parse_args()
	print args

	t=isodelta.database()
	t.read(args.database)
	
if __name__=="__main__":
	main()
