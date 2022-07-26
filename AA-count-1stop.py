#!/usr/local/bin/python

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy


# set some definitions

seqLen = 11
aaList = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*"]

# Import data

records = list(SeqIO.parse("SampleName.translate.fasta", "fasta"))

# For each sequence with 11 residues, iterate through the positions, determine the amino acid identity, and add a count to the appropriate position on the frequency matrix.

frequencyMatrix = numpy.zeros(shape=(len(aaList),seqLen))
for i in range(len(records)):
	if len(records[i].seq) == seqLen and records[i].seq[5] == "Y" and records[i].seq.count("*") == 1:
		for pos in range(seqLen):
			for aa in range(len(aaList)):
				if records[i].seq[pos] == aaList[aa]:
					frequencyMatrix[aa,pos] +=1

# Save the array as a tab-delimited text file.

numpy.savetxt("SampleName.AA-count-1stop.txt", frequencyMatrix, fmt = '%i', delimiter='\t')
