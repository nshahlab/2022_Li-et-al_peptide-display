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
		peptideList = list()
		peptideListVar = OrderedDict()
		peptideListRef = OrderedDict()

		with open("pTyr-Var_variant.txt") as f:
			for line in f:
				peptide = line.strip('\n\r')
				peptideList.append([peptide,0,str(),0])
				peptideListVar[peptide] = 0

		with open("pTyr-Var_reference.txt") as f:
			counter = 0
			for line in f:
				peptide = line.strip('\n\r')
				peptideList[counter][2] = peptide
				counter += 1
				peptideListRef[peptide] = 0

		for nuc in SeqIO.parse(filename, "fastq"):
			translation = SeqRecord(seq = nuc.seq.translate(to_stop=False))
			if str(translation.seq) in peptideListVar:
				peptideListVar[str(translation.seq)] += 1
			if str(translation.seq) in peptideListRef:
				peptideListRef[str(translation.seq)] += 1

		for pair in peptideList:
			if pair[0] in peptideListVar:
				pair[1] = peptideListVar[pair[0]]
			if pair[2] in peptideListRef:
                		pair[3] = peptideListRef[pair[2]]

		countsFile = open(filename[0:-13] + 'peptide_counts-var-ref.txt', 'w+')
		for pair in peptideList:
			countsFile.write(pair[0] + "\t" + str(pair[1]) + "\t" + pair[2] + "\t" + str(pair[3]) + "\n")

		countsFile.close()
