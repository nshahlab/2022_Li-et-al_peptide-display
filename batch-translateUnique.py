#!/usr/local/bin/python

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter
import subprocess
import os
from collections import OrderedDict

# define a variable that is the string for the present working directory

directory = subprocess.check_output(['pwd'], universal_newlines=True).rstrip()

# This first part will read a trimmed fastq file of your coding sequence and output a translated fasta file

def make_protein_record(nuc_record):
	return SeqRecord(seq = nuc_record.seq.translate(to_stop=False), \
		id = "trans_" + nuc_record.id, \
			description = "translation of MiSeq read")

for filename in os.listdir(directory):
	if filename.endswith(".trimmed.fastq"):
		translatedFile = filename[0:-13] + "translate.fasta"	
		proteins = (make_protein_record(nuc_rec) for nuc_rec in \
			SeqIO.parse(filename, "fastq"))
		writer = FastaWriter(open(translatedFile, "w"), wrap=0)
		writer.write_file(proteins)
		proteins = None

# This part will make an ordered dictionary listing every unique peptide of length = seqLen encoded in the trimmed fastq file and count occurrences (note that the first part of this script isn't needed for the second part).

seqLen = 11
stop = "*"

for filename in os.listdir(directory):
	if filename.endswith(".translate.fasta"):
		peptideList = OrderedDict()
		nonStopList = OrderedDict()
		stopList = OrderedDict()
		for translation in SeqIO.parse(filename, "fasta"):
			if not str(translation.seq) in peptideList and len(str(translation.seq)) == seqLen:
				peptideList[str(translation.seq)] = 1
			elif len(str(translation.seq)) == seqLen:
				peptideList[str(translation.seq)] += 1
		for peptide, count in peptideList.items():
			residueList = list(peptide)
			if not stop in residueList:
				nonStopList[peptide] = count
			elif stop in residueList:
				stopList[peptide] = count
		fullCounts = open(filename[0:-15] + 'peptide-counts_full.txt', 'w+')
		for peptide in peptideList:
			fullCounts.write(peptide + "\t" + str(peptideList[peptide]) + "\n")
		fullCounts.close()
		nonStopCounts = open(filename[0:-15] + 'peptide-counts_non-stop.txt', 'w+')
		for peptide in nonStopList:
			nonStopCounts.write(peptide + "\t" + str(peptideList[peptide]) + "\n")
		nonStopCounts.close()
		stopCounts = open(filename[0:-15] + 'peptide-counts_stop.txt', 'w+')
		for peptide in stopList:
			stopCounts.write(peptide + "\t" + str(peptideList[peptide]) + "\n")
		stopCounts.close()
