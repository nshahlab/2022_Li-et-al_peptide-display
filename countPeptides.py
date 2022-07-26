#!/usr/local/bin/python

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import numpy
from collections import OrderedDict

# Import peptide list

peptideList = OrderedDict()
with open("pTyr-Var_full.txt") as f:
    for line in f:
        peptide = line.strip('\n\r')
        peptideList[peptide] = 0

# Count peptide frequency

for nuc in SeqIO.parse("SampleName.trimmed.fastq", "fastq"):
    translation = SeqRecord(seq = nuc.seq.translate(to_stop=False))
    if str(translation.seq) in peptideList:
        peptideList[str(translation.seq)] += 1        

countsFile = open('SampleName.peptide-counts.txt', 'w+')
for peptide in peptideList:
    countsFile.write(peptide + "\t" + str(peptideList[peptide]) + "\n")
countsFile.close()
