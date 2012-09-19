#!/usr/bin/python

import isodelta
import argparse

# Update database for energy shifts in EuO and substrate
def main():
	parser = argparse.ArgumentParser(description='Update database for energy shifts in EuO and substrate')
	parser.add_argument('folders', type=str, nargs='+', help='Folders containing results of isolated material runs')
  	parser.add_argument('-d', '--database', type=str, help='Database file name', default='isodelta.dat')

	args = parser.parse_args()

	t=isodelta.database()
	t.extract(args.folders)
	t.write(args.database)
	
if __name__=="__main__":
	main()
