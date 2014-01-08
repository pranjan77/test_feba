
# Jim Kaminski
# Arkin Lab


import sys
import csv
import argparse
from argparse import RawTextHelpFormatter
import subprocess
import re
import os
import time
import datetime
import math

parser = argparse.ArgumentParser(description='This program counts up the insertions in each gene.')


grpInput = parser.add_argument_group('Input')
grpInput.add_argument('--genes', type=str, dest='strGenes',default= "", help='Enter the list of genes.')
grpInput.add_argument('--pool', type=str, dest='strInserts',default= "", help='Enter the pool information.')
args = parser.parse_args()
# Genes can be on the megaplasmid or the chromosome. Make a dictionary for each.
dictChromosome = {}
dictMegaplasmid = {}




"""
Genes.tab
Field, Index
locusId 0
sysName 1
type    2
scaffoldId  3
begin   4
end 5
strand  6
name    7
desc    8
"""
#################################################################################
# Load in the genes.tab file. Make two sorted lists of positions. One for the
# Chromosome and one for the Megaplasmid.
#################################################################################
with open(args.strGenes, 'rb') as fileGenes:
	fileGenes.next()
	aaChromosome = []
	aaMegaplasmid = []
	for astrLine in csv.reader(fileGenes,delimiter='\t'):
		iScaffoldId = int(astrLine[3])
		strSysName = astrLine[1]
		iStart = int(astrLine[4])
		iEnd = int(astrLine[5])

		if iScaffoldId == 139:
			aaChromosome.append([iStart,iEnd,0,0] + astrLine)
		elif iScaffoldId == 140:
			aaMegaplasmid.append([iStart,iEnd,0,0] + astrLine)
		else:
			sys.stderr.write("Error: Scaffold ID= " + str(iScaffoldId) + " Cannot determine if chromosome or megaplasmid.\n")

# Sort the lists

aaChromosome.sort(key=lambda x: float(x[0]))
aaMegaplasmid.sort(key=lambda x: float(x[0]))

################################################################################
# Load in the Insertions
################################################################################

"""
barcode 0
rcbarcode   1
nTot    2
n   3
scaffold    4
strand  5
pos 6
n2  7
scaffold2   8
strand2 9
pos2    10
nPastEnd    11

"""

################################################################################
# Binary Search
# What if the end part lands in the region of interest? (iPos + 19?)
# What i

# -/+ - orientation of the transposon.
# ignore

################################################################################
def FindClosestPos(aaChromosome, iPos):
	# Check the middle value
	iStart = 0
	iEnd = len(aaChromosome)-1
	iIndex = (iStart+iEnd)/2

	while (iStart != iEnd):
		#print iStart, iEnd, iIndex, aaChromosome[iIndex][0], aaChromosome[iIndex][1]

		if (aaChromosome[iIndex][0] <= iPos and aaChromosome[iIndex][1] >= iPos) or (aaChromosome[iIndex][0] <= (iPos + 19) and aaChromosome[iIndex][1] >= (iPos + 19)) :
			return iIndex

		elif aaChromosome[iIndex][0] > (iPos):
			iStart = iStart
			iEnd = iIndex
			iIndex = (iStart + iEnd)/2

		else:
			iStart = iIndex +1
			iEnd = iEnd
			iIndex = (iStart + iEnd)/2

	if aaChromosome[iIndex][0] <= iPos and aaChromosome[iIndex][1] >= iPos:
		return iIndex
	else:
		return -1


with open(args.strInserts, 'rb') as fileInserts:
	fileInserts.next()
	iGene = 0
	for astrLine in csv.reader(fileInserts,delimiter='\t'):

		aGeneOne = astrLine[0:3] + astrLine[3:7] + [astrLine[11]]
		aGeneTwo = astrLine[0:3] + astrLine[7:11] + [astrLine[11]]
		# aGeneOne = [barcode,rcbarcode,nTot,n,scaffold,strand,pos,nPastEnd]
        # aGeneTwo = [barcode,rcbarcode,nTot,n2,scaffold2,strand2,pos2,nPastEnd]

		# We ignore aGeneTwo for now. That was just a way of checking instances
		# where we found the transposon elsewhere, which shouldn't happen.
		# This was just due to chimeric PCR.

		for aGene in [aGeneOne]:
			if aGene[4] != "pastEnd" and aGene[4] != "":
				#print aGene[3:7]
				iGene+=1
				iScaffoldId = int(aGene[4])
				if aGene[5]=="+":
					iStrand = 2
				elif aGene[5]=="-":
					iStrand = 3

				iPos = int(aGene[6])

				#print iScaffoldId, iPos, aGene[5]

				if iScaffoldId==139:
					iIndex = FindClosestPos(aaChromosome, iPos)
					if iIndex>=0:
						aaChromosome[iIndex][iStrand] = aaChromosome[iIndex][iStrand]+int(aGene[3])
				if iScaffoldId==140:
					iIndex = FindClosestPos(aaMegaplasmid, iPos)
					if iIndex>=0:
						aaMegaplasmid[iIndex][iStrand] = aaMegaplasmid[iIndex][2]+int(aGene[3])


print "Chromosome"
for aGene in aaChromosome:
	#iTotal = aGene[1] +aGene[2]
	print '\t'.join(map(str,aGene))

print "Megaplasmid"
for aGene in aaMegaplasmid:
	print '\t'.join(map(str,aGene))
