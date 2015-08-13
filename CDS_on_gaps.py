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


def find_CDS_on_Gaps(fileName):

	prevName = ''
	firstTime = True
	intervalToCheck = [-1,-1, '']
	CDS_in_gap = []

	with open(fileName, 'r') as gff:
		for i in gff:
			if '#' not in i:
				splitLine = i.split('\t')
				typeRegion = splitLine[2]
				name = splitLine[0]
				if firstTime:
					prevName = name
					intervalToCheck[2] = name

				if 'gap' in typeRegion:
					intervalToCheck[0] = int(splitLine[3])
					intervalToCheck[1] = int(splitLine[4])
					if prevName != name:
						intervalToCheck[2] = name
					prevName = name

				if 'CDS' in typeRegion and name == intervalToCheck[2] and intervalToCheck[0] <= int(splitLine[3]) and intervalToCheck[1] >= int(splitLine[4]):
					CDS_in_gap.append(i)

	return CDS_in_gap


def write_CDS_in_gap_file(CDSarray, fileName):
	
	with open(fileName+'_CDS_in_gaps.gff', 'w') as gapFile:
		for i in CDSarray:
			gapFile.write(i)


