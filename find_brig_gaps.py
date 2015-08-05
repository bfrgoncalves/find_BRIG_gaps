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


def main():

	parser = argparse.ArgumentParser(description="This program finds the gaps of a query genome after the alignment againsta reference using BRIG software.")
	parser.add_argument('-x', nargs='?', type=str, help="query .tab file", required=True)
	parser.add_argument('-ib', nargs='?', type=str, help="interval begin", required=True)
	parser.add_argument('-ie', nargs='?', type=str, help="interval end", required=True)
	parser.add_argument('-s', nargs='?', type=str, help="sensitivity", required=True)
	parser.add_argument('-f', nargs='?', type=str, help="fasta from sequences used as reference on BRIG", required=True)
	parser.add_argument('-o', nargs='?', type=str, help="results folder", required=True)


	args = parser.parse_args()


	LineDict, arrayOfLines = importTsv(args.x, '\t')

	orderedTab = orderTsv(LineDict, arrayOfLines)

	gaps = getGaps(orderedTab, args.ib, args.ie, args.s)

	gapArray = getGapSequence(args.f, gaps)
	writeGapFile(args.o, gapArray)

	print gaps


def importTsv(fileName, delimiter):

	LineDict = {}
	arrayOfLines = []

	with open(fileName) as tsvFile:
		for line in csv.reader(tsvFile, delimiter = delimiter): #You can also use delimiter="\t" rather than giving a dialect.

			if float(line[9]) < float(line[8]):
				begin = line[9]
				end = line[8]
				
			else:
				begin = line[8]
				end = line[9]

			identifier = begin + '-' + end

			if identifier not in arrayOfLines:
				arrayOfLines.append((identifier, float(begin), float(end)))
				LineDict[identifier] = line


	return LineDict, arrayOfLines


def orderTsv(LineDict, arrayOfLines):

	arrayOfLines = sorted(arrayOfLines, key=lambda tuple: tuple[1])
	
	return arrayOfLines


def getGaps(orderedTab, Ibegin, Iend, sense):

	gaps = []
	coveredRegion = 0

	for i in range(0, len(orderedTab)-2):
		gap = orderedTab[i][1] - coveredRegion
		#print gap
		if float(Ibegin) <= orderedTab[i][1] and float(Iend) >= orderedTab[i][1]:
			if gap >= float(sense):
				gaps.append((str(coveredRegion) + '--' + str(coveredRegion + gap), gap))

		if gap > 0:
			if gap + orderedTab[i][2] > coveredRegion:
				coveredRegion += (gap + orderedTab[i][2] - coveredRegion)
		elif orderedTab[i][2] - coveredRegion > 0:
			coveredRegion += orderedTab[i][2] - coveredRegion

	print coveredRegion
	return gaps


def getGapSequence(fastaFile, gaps):

	gene_fp = HTSeq.FastaReader(fastaFile)
	wholeSequence = ''
	newGapFile = []
	for allele in gene_fp:
		wholeSequence += allele.seq

	for i in gaps:
		begin = float(i[0].split('--')[0])
		end = float(i[0].split('--')[1])
		seqToExtract = wholeSequence[int(begin):int(end)]
		newGapFile.append((i[0],i[1],seqToExtract))

	return newGapFile

def writeGapFile(fileName, gapArray):

	with open(fileName, 'w') as gapFile:
		for i in gapArray:
			gapFile.write('>region_' + str(i[0]) + '_gapLength_' + str(int(float(i[1]))) + '\n' + i[2] + '\n')



if __name__ == "__main__":
	main()