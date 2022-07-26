#!/usr/local/bin/python

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict
import subprocess
import os

directory = subprocess.check_output(['pwd'], universal_newlines=True).rstrip()

for filename in os.listdir(directory):
	if filename.endswith(".trimmed.fastq"):
		peptideList = OrderedDict()
		with open("pTyr-Var_full.txt") as f:
			for line in f:
				peptide = line.strip('\n\r')
				peptideList[peptide] = 0

		for nuc in SeqIO.parse(filename, "fastq"):
			translation = SeqRecord(seq = nuc.seq.translate(to_stop=False))
			if str(translation.seq) in peptideList:
				peptideList[str(translation.seq)] += 1        

		countsFile = open(filename[0:-13] + 'peptide_counts.txt', 'w+')
		for peptide in peptideList:
			countsFile.write(peptide + "\t" + str(peptideList[peptide]) + "\n")
		countsFile.close()
