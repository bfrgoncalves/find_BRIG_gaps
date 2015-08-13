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

from BRIG_gaps_utils import runPipeline


def main():

	parser = argparse.ArgumentParser(description="This program finds the gaps of a query genome after the alignment against a reference using BRIG software.")
	parser.add_argument('-x', nargs='?', type=str, help="query .tab file from BRIG scratch folder/folder with .tab scratch files", required=True)
	parser.add_argument('-ib', nargs='?', type=str, help="interval begin", required=True)
	parser.add_argument('-ie', nargs='?', type=str, help="interval end", required=True)
	parser.add_argument('-s', nargs='?', type=str, help="sensitivity", required=True)
	parser.add_argument('-n', nargs='?', type=str, help="Gaps name", required=True)
	parser.add_argument('-m', nargs='?', type=str, help="Is multiple Contigs file", required=True)
	parser.add_argument('-f', nargs='?', type=str, help="fasta from sequences used as reference on BRIG", required=True)
	parser.add_argument('-a', nargs='?', type=str, help=".gff file to change", required=False)
	parser.add_argument('-o', nargs='?', type=str, help="results folder", required=True)


	args = parser.parse_args()

	if os.path.isdir(args.x):
		onlyfiles = [ f for f in listdir(args.x) if isfile(join(args.x,f)) ]
		if args.a:
			gffFile = args.a
			gffBasename = os.path.basename(gffFile)
			gffName = os.path.splitext(gffBasename)
		for i in onlyfiles:
			gapName = os.path.splitext(i)[0] + '_Gap'
			resultsFileName = args.o + '_' + os.path.splitext(i)[0]
			runPipeline(os.path.join(args.x,i), args.ib, args.ie, args.s, gapName, args.f, gffFile, resultsFileName, True, args.m)

		if args.a:
			os.rename(gffName[0]+'_inter.gff', gffName[0]+'_plusGAPS.gff')
	else:
		runPipeline(args.x, args.ib, args.ie, args.s, args.n, args.f, args.a, args.o, False, args.m)


if __name__ == "__main__":
	main()