import subprocess
import argparse
import os
import shutil
from os import listdir
from os.path import isfile, join, isdir
import sys
from datetime import datetime
import numpy
import string
import operator
import collections
import csv
import HTSeq

from Bio.Blast import NCBIXML


def main():

	parser = argparse.ArgumentParser(description="This program parses a .xml BLAST results file into a .tab BRIG format file")
	parser.add_argument('-x', nargs='?', type=str, help=".xml BLAST results file", required=True)
	parser.add_argument('-o', nargs='?', type=str, help="results file name", required=True)


	args = parser.parse_args()

	BLAST_to_BRIG(args.x, args.o)


def BLAST_to_BRIG(BLASTfile, resultsFile):

	rec = open(BLASTfile)
	blast_records = NCBIXML.parse(rec)

	with open(resultsFile, 'w') as tabFile: 
	
		for blast_record in blast_records:
		
			for alignment in blast_record.alignments:
				for match in alignment.hsps:
					tabFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
								blast_record.query,
								alignment.hit_def,
								round(float(match.identities)/float(alignment.length),2),
								int(match.score),
								alignment.length,
								int(alignment.length) - int(match.identities),
								match.query_start,
								(int(match.query_start) + int(alignment.length)),
								match.sbjct_start,
								(int(match.query_start) + int(alignment.length))))

					break


if __name__ == "__main__":
	main()