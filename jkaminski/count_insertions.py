###############################################################################
# Jim Kaminski
# Arkin Lab
# Dec 2013 - Jan 2014
###############################################################################

import sys
import csv
import argparse
from argparse import RawTextHelpFormatter
import os
import logging
import datetime
import math

logging.basicConfig(filename='genes.log',level=logging.DEBUG)
logging.info(datetime.date.today().ctime())

parser = argparse.ArgumentParser(description='This program counts up the insertions in each gene.')

grpInput = parser.add_argument_group('Input')
grpInput.add_argument('--genes', type=str, dest='strGenes',default= "", help='Enter the list of genes.')
grpInput.add_argument('--pool', type=str, dest='strInserts',default= "", help='Enter the pool information.')
grpInput.add_argument('--pctoff', type=float, dest='dOffset',default= 0.0, help='This will ignore insertions in the first and last x% of the gene (.10,.20, etc.) ')
args = parser.parse_args()



"""
Genes.tab - Header Information:
(Field,Index)

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
# Load in the genes.tab file. For each unique ScaffoldId, make a sorted list of
# positions, and save these in dictScaffoldPosLists.
#################################################################################

dictScaffoldPosLists = {}

with open(args.strGenes, 'rb') as fileGenes:
	fileGenes.next()

	for astrLine in csv.reader(fileGenes,delimiter='\t'):
		strScaffoldID = str(astrLine[3])
		strSysName = astrLine[1]
		iStart = int(astrLine[4])
		iEnd = int(astrLine[5])
		strDesc = astrLine[8]

		if strScaffoldID not in dictScaffoldPosLists:
			dictScaffoldPosLists[strScaffoldID] = []

		dictScaffoldPosLists[strScaffoldID] = dictScaffoldPosLists[strScaffoldID] + [[iStart,iEnd,0,0]+[strSysName]+[strDesc]]

# Sort the lists
for key in dictScaffoldPosLists.keys():
	dictScaffoldPosLists[key] = sorted(dictScaffoldPosLists[key], key=lambda x: float(x[0]))


# DEBUGGING - Test to ensure the dictionaries are correct.
"""
for key in dictScaffoldPosLists.keys():
	print key
	for x in dictScaffoldPosLists[key]:
		print x
"""

################################################################################
# Binary Search - Check if Insertion landed in a gene.

# Look into later:
#   Genes overlap, and I currently just get the first hit. To solve this problem
#       Return an array of hits.
#       If you get a hit, check next gene below it. Add to list of hits if correct. If not stop.
#       Then check next gene above hit. Add to list of hits if it's correct. If not, stop.

#   What if the end part lands in the region of interest? (iPos + 19?)
################################################################################
def FindClosestPos(aaChromosome, iPos,dPctOffset):
	# Check the middle value
	iStart = 0
	iEnd = len(aaChromosome)-1
	iIndex = (iStart+iEnd)/2
	aiHits = []


	while (iStart != iEnd and not aiHits):

		iOffset = math.floor(dPctOffset * (aaChromosome[iIndex][1] - aaChromosome[iIndex][0]))
		iStartGene = aaChromosome[iIndex][0] + iOffset
		iEndGene = aaChromosome[iIndex][1] - iOffset

		#print iOffset, iStartGene, iEndGene

		if (iStartGene <= iPos and iEndGene >= iPos) or (iStartGene <= (iPos + 19) and iEndGene >= (iPos + 19)) :
			aiHits = aiHits + [iIndex]

		elif iStartGene > (iPos):
			iStart = iStart
			iEnd = iIndex
			iIndex = (iStart + iEnd)/2

		else:
			iStart = iIndex +1
			iEnd = iEnd
			iIndex = (iStart + iEnd)/2

	if not aiHits:
		return [-1]


	# If you have a hit, check below and above for more hits.
	bCheckedLow = False
	bCheckedHigh = False

	iLow = iIndex -1
	iHigh = iIndex + 1

	# Check the genes below the hit.
	while (bCheckedLow == False and iLow>0):

		iStartLow = aaChromosome[iLow][0] + iOffset
		iEndLow = aaChromosome[iLow][1] - iOffset

		if iStartLow <= iPos and iEndLow >= iPos:
			aiHits = aiHits + [iLow]
		else:
			bCheckedLow = True
		iLow = iLow-1

	while (bCheckedHigh == False and iHigh<len(aaChromosome)):

		iStartHigh = aaChromosome[iHigh][0] + iOffset
		iEndHigh = aaChromosome[iHigh][1] - iOffset

		if aaChromosome[iHigh][0] <= iPos and aaChromosome[iHigh][1] >= iPos:
			aiHits = aiHits + [iHigh]
		else:
			bCheckedHigh = True
		iHigh = iHigh+1

	return aiHits





"""
Pool.n10 - Header Information:
(Field,Index)

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
# Load in the Insertions
################################################################################
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
				strScaffoldId = str(aGene[4])
				if aGene[5]=="+":
					iStrand = 2
				elif aGene[5]=="-":
					iStrand = 3

				iPos = int(aGene[6])

				aiIndex = FindClosestPos(dictScaffoldPosLists[strScaffoldId], iPos, args.dOffset)
				for iIndex in aiIndex:
					if iIndex>=0:
						dictScaffoldPosLists[strScaffoldId][iIndex][iStrand] = dictScaffoldPosLists[strScaffoldId][iIndex][iStrand]+int(aGene[3])

iNoInsertions = 0

for key in dictScaffoldPosLists.keys():
	for aGene in dictScaffoldPosLists[key]:
		iCount = aGene[2] + aGene[3]
		print '\t'.join([aGene[4],str(iCount),aGene[-1]])
		#print '\t'.join(map(str,aGene))
		if iCount==0:
			iNoInsertions+=1

print "Genes_without_insertions:\t"+str( iNoInsertions)
