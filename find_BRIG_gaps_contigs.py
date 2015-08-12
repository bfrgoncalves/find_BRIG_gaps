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

	parser = argparse.ArgumentParser(description="This program finds the gaps of a query genome after the alignment against a reference using BRIG software.")
	parser.add_argument('-x', nargs='?', type=str, help="query .tab file from BRIG scratch folder/folder with .tab scratch files", required=True)
	parser.add_argument('-ib', nargs='?', type=str, help="interval begin", required=True)
	parser.add_argument('-ie', nargs='?', type=str, help="interval end", required=True)
	parser.add_argument('-s', nargs='?', type=str, help="sensitivity", required=True)
	parser.add_argument('-n', nargs='?', type=str, help="Gaps name", required=True)
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
			runPipeline(os.path.join(args.x,i), args.ib, args.ie, args.s, gapName, args.f, gffFile, resultsFileName, True)

		if args.a:
			os.rename(gffName[0]+'_inter.gff', gffName[0]+'_plusGAPS.gff')
	else:
		runPipeline(args.x, args.ib, args.ie, args.s, args.n, args.f, args.a, args.o, False)


def runPipeline(x, ib, ie, s, n, f, a, o, overwrite):

	LineDict, arrayOfLines = importTsv(x, '\t')
	print
	orderedTab = orderTsv(arrayOfLines)
	#print orderedTab
	gaps, totalCoveredZone = getGaps(orderedTab, ib, ie, s)
	gapArray, wholeSequenceSize = getGapSequence(f, gaps)
	#print gapArray
	print 'Coverage percentage: ' + str((float(totalCoveredZone)/float(wholeSequenceSize)) * 100) + '%'
	print 
	
	if n:
		gapName = n
	else:
		gapName = 'Gap'
	
	writeGapFiles(o, gapArray, f, gapName)

	if a:
		changeGff(a, gapArray, gapName, overwrite)

	print 'Gaps: ' + str(gaps) + '\n'



def importTsv(fileName, delimiter):

	LineDict = {}
	arrayOfLines = []

	with open(fileName) as tsvFile:
		for line in csv.reader(tsvFile, delimiter = delimiter): #You can also use delimiter="\t" rather than giving a dialect.

			if float(line[9]) < float(line[8]):
				contig_Align_begin = line[9]
				contig_Align_end = line[8]
				
			else:
				contig_Align_begin = line[8]
				contig_Align_end = line[9]

			referenceSplit = line[1].split(':')
			reference = referenceSplit[0]
			referenceBeginPlace = referenceSplit[1]
			referenceEndPlace = referenceSplit[2]

			identifier = referenceBeginPlace + '_' + referenceEndPlace + '_' + contig_Align_begin + '_' + contig_Align_end

			if identifier not in arrayOfLines:
				arrayOfLines.append((reference, int(referenceBeginPlace), int(referenceEndPlace), identifier, int(contig_Align_begin), int(contig_Align_end)))
				LineDict[identifier] = line


	return LineDict, arrayOfLines


def orderTsv(arrayOfLines):
	newArrayOfLines = []
	intermediaryArray = []

	arrayOfLines = sorted(arrayOfLines, key=lambda tuple: tuple[1])

	prevReference = ''
	for i in arrayOfLines:
		if prevReference == i[0]:
			intermediaryArray.append(i)
		else:
			intermediaryArray = sorted(intermediaryArray, key=lambda tuple: tuple[4])
			for j in intermediaryArray:
				newArrayOfLines.append(j)
			intermediaryArray = []
			prevReference = i[0]
			intermediaryArray.append(i)
	
	return newArrayOfLines


def getGaps(orderedTab, Ibegin, Iend, sense):

	gaps = []
	coveredRegion = 0
	totalCoveredZone = 0
	totalGapSize = 0
	prevRefName = orderedTab[0][0]
	prevContigEnd = orderedTab[0][2]
	BeginOfContig = 0

	for i in range(0, len(orderedTab)-2):
		if prevRefName != orderedTab[i][0]:
			BeginOfContig = orderedTab[i][1]
			prevRefName = orderedTab[i][0]
			coveredRegion = prevContigEnd
		else:
			prevContigEnd = orderedTab[i][2]
	
		gap = BeginOfContig + orderedTab[i][4] - coveredRegion
		#print gap
		#print BeginOfContig, coveredRegion, gap
		if float(Ibegin) <= coveredRegion and float(Iend) >= coveredRegion + gap:
			if gap >= float(sense):
				gaps.append((str(coveredRegion-BeginOfContig) + '--' + str(coveredRegion - BeginOfContig + gap), gap, orderedTab[i][0]))

		if gap > 0:
			totalGapSize += gap
			if gap + BeginOfContig + orderedTab[i][5] > coveredRegion:
				coveredRegion += (gap + (BeginOfContig + orderedTab[i][5] - (gap + coveredRegion)))
		elif BeginOfContig + orderedTab[i][5] - coveredRegion > 0:
			coveredRegion += BeginOfContig + orderedTab[i][5] - coveredRegion
	
	totalCoveredZone += coveredRegion - totalGapSize

	print 'Covered Region: ' + str(int(totalCoveredZone)) + ' bp'
	return gaps, totalCoveredZone


def getGapSequence(fastaFile, gaps):

	gene_fp = HTSeq.FastaReader(fastaFile)
	wholeSequence = {}
	newGapFile = []
	for allele in gene_fp:
		try:
			wholeSequence[allele.name] += allele.seq
		except KeyError:
			wholeSequence[allele.name] = allele.seq

	print 'Reference Size: ' + str(len(wholeSequence)) + ' bp'

	for i in gaps:
		begin = float(i[0].split('--')[0])
		end = float(i[0].split('--')[1])
		seqToExtract = wholeSequence[i[2]][int(begin):int(end)]
		newGapFile.append((i[0],i[1],seqToExtract, i[2]))

	return newGapFile, len(wholeSequence)

def writeGapFiles(fileName, gapArray, fastaName, gapName):

	with open(fileName+'.fasta', 'w') as gapFile:
		for i in gapArray:
			gapFile.write('>region_' + str(i[3]) + '_' + str(i[0]) + '_gapLength_' + str(int(float(i[1]))) + '\n' + i[2] + '\n')


	with open(fileName+'.gff', 'w') as gapFile:
		gapFile.write('##gff\n##dataFROM_find_brig_gaps_@bfrgoncalves\n')
		fastaName = os.path.basename(fastaName)
		countGaps = 0
		for i in gapArray:
			countGaps += 1
			gapFile.write(i[3] + '\tRefSeq\tgap\t' + str(int(float(i[0].split('--')[0]))) +'\t' + str(int(float(i[0].split('--')[1]))) + '\t.\t.\t.\t' + 'ID=' + gapName + str(countGaps) + ';Name=' + gapName + str(countGaps) + '_' + str(int(i[1])) + '\n')


def changeGff(gffFile, gapArray, gapName, overwrite):
	
	gffBasename = os.path.basename(gffFile)
	gffName = os.path.splitext(gffBasename)
	currentGapArray = gapArray
	countGaps = 0
	firstTime = True

	if os.path.isfile(gffName[0]+'_inter.gff'):
		gffFile = gffName[0] + '_inter.gff'
	
	with open(gffFile, 'r') as gff:
		with open(gffName[0]+'_plusGAPS.gff', 'w') as gapFile:
			for i in gff:
				if '#' not in i:
					if firstTime:
						sequenceName = i.split('\t')[0]
						firstTime = False
					try:
						gffBegin = float(i.split('\t')[3])
						referenceName = i.split('\t')[0]
					except IndexError:
						gapFile.write(i)
						continue
					if not currentGapArray:
						gapFile.write(i)
					else:
						currentGap = currentGapArray[0]
						nextGapBegin = float(currentGap[0].split('--')[0])

						if nextGapBegin > gffBegin or referenceName != currentGap[3]:
							gapFile.write(i)
						else:
							countGaps += 1
							gapFile.write(currentGap[3] + '\tRefSeq\tgap\t' + str(int(float(currentGap[0].split('--')[0]))) +'\t' + str(int(float(currentGap[0].split('--')[1]))) + '\t.\t.\t.\t' + 'ID=' + gapName + str(countGaps) + ';Name=' + gapName + str(countGaps) + '_' + str(int(currentGap[1])) + '\n')
							currentGapArray.pop(0)

				else:
					gapFile.write(i)

	if overwrite:
		if os.path.isfile(gffName[0]+'_inter.gff'):
			os.remove(gffName[0]+'_inter.gff')
			os.rename(gffName[0]+'_plusGAPS.gff', gffName[0]+'_inter.gff')
		else:
			os.rename(gffName[0]+'_plusGAPS.gff', gffName[0]+'_inter.gff')


if __name__ == "__main__":
	main()