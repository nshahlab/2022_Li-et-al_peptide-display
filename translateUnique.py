#!/usr/local/bin/python

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter
from collections import OrderedDict

# This first part will read a trimmed fastq file of your coding sequence and output a translated fasta file

def make_protein_record(nuc_record):
	return SeqRecord(seq = nuc_record.seq.translate(to_stop=False), \
		id = "trans_" + nuc_record.id, \
			description = "translation of MiSeq read")

from Bio import SeqIO
proteins = (make_protein_record(nuc_rec) for nuc_rec in \
	SeqIO.parse("SampleName.trimmed.fastq", "fastq"))
writer = FastaWriter(open("SampleName.translate.fasta", "w"), wrap=0)
writer.write_file(proteins)

# This part will make an ordered dictionary listing every unique peptide of length = seqLen encoded in the trimmed fastq file and count occurrences (note that the first part of this script isn't needed for the second part).

seqLen = 11
peptideList = OrderedDict()

for nuc in SeqIO.parse("SampleName.trimmed.fastq", "fastq"):
	translation = SeqRecord(seq = nuc.seq.translate(to_stop=False))
	if not str(translation.seq) in peptideList and len(str(translation.seq)) == seqLen:
		peptideList[str(translation.seq)] = 0
		peptideList[str(translation.seq)] += 1
	elif len(str(translation.seq)) == seqLen:
		peptideList[str(translation.seq)] += 1

# this part creates similar order dictionaries but that are separated based on whether or not the translated sequences of length seqLen have a stop codon or not

nonStopList = OrderedDict()
stopList = OrderedDict()
stop = "*"

for peptide, count in peptideList.items():
	residueList = list(peptide)
	if not stop in residueList:
		nonStopList[peptide] = count
	elif stop in residueList:
		stopList[peptide] = count

# everything is written to files

fullCounts = open('SampleName.peptide-counts_full.txt', 'w+')
for peptide in peptideList:
    fullCounts.write(peptide + "\t" + str(peptideList[peptide]) + "\n")
fullCounts.close()

nonStopCounts = open('SampleName.peptide-counts_non-stop.txt', 'w+')
for peptide in nonStopList:
    nonStopCounts.write(peptide + "\t" + str(peptideList[peptide]) + "\n")
nonStopCounts.close()

stopCounts = open('SampleName.peptide-counts_stop.txt', 'w+')
for peptide in stopList:
    stopCounts.write(peptide + "\t" + str(peptideList[peptide]) + "\n")
stopCounts.close()

