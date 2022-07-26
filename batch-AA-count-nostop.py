#!/usr/local/bin/python

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy
import subprocess
import os

# define a variable that is the string for the present working directory

directory = subprocess.check_output(['pwd'], universal_newlines=True).rstrip()

# set some definitions

seqLen = 11
aaList = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*"]

# Import data

for filename in os.listdir(directory):
	if filename.endswith(".translate.fasta"):
		records = list(SeqIO.parse(filename, "fasta"))
		frequencyMatrix = numpy.zeros(shape=(len(aaList),seqLen))
		for i in range(len(records)):
			if len(records[i].seq) == seqLen and records[i].seq[5] == "Y" and "*" not in records[i].seq:
				for pos in range(seqLen):
					for aa in range(len(aaList)):
						if records[i].seq[pos] == aaList[aa]:
							frequencyMatrix[aa,pos] +=1
		numpy.savetxt(filename[0:-15] + 'AA-count-nostop.txt', frequencyMatrix, fmt = '%i', delimiter='\t')
