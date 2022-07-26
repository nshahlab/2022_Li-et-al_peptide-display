#Inputs: 
#Counts table for whatever selected sample you are interested in. 
#Counts table for input library. 
#List of peptides to score in a file called peptides.txt

import numpy

# Define sequence length and amino acid list order.

seqLen = 11
aaList = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*"]

#specify the counts tables you are interested in, below.

dataTest = numpy.loadtxt("SelectedSampleName.AA-count-nostop.txt", delimiter='\t')
dataRef = numpy.loadtxt("InputSampleName.AA-count-nostop.txt", delimiter='\t')
dataTestFreq = dataTest[0:20]/(dataTest.sum()/seqLen)
dataRefFreq = dataRef[0:20]/(dataRef.sum()/seqLen)
dataNorm = (dataTestFreq/dataRefFreq)

#this section will determine the minimum and maximum score possible, to be used for normalization

maxScore = 0
minScore = 0

for i in range(seqLen):
    if i !=5:
        maxScore = maxScore + numpy.log2(numpy.max(dataNorm[:,i]))
        minScore = minScore + numpy.log2(numpy.min(dataNorm[:,i]))

maxScoreFinal = maxScore/(seqLen-1)
minScoreFinal = minScore/(seqLen-1)
print("max: " + str(maxScoreFinal))
print("min: " + str(minScoreFinal))

#this section reads the input list of sequences and also creates two output files.

with open('peptides.txt') as f:
    peptides = f.readlines()

scoreFile = open('SelectedSampleName_peptide-scores.txt', 'w+')
scoreFileNorm = open('SelectedSampleName_peptide-scores-norm.txt', 'w+')
 
# calculates the score for each peptide in the input file. The scoring is only done for positions -5 to -1 and 1 to 5. 

for pep in peptides: 
    score=0
    for i in range(len(pep)):
        if i !=5:
            for aa in range(len(aaList)):
                if pep[i]==aaList[aa]:                
                    score= score+ numpy.log2(dataNorm[aa,i])
    scoreFinal = score/(seqLen-1)
    scoreFinalNorm = (scoreFinal - minScoreFinal)/(maxScoreFinal - minScoreFinal)
    scoreFile.write(pep[0:11] + "\t" + str(scoreFinal) + "\n")
    scoreFileNorm.write(pep[0:11] + "\t" + str(scoreFinalNorm) + "\n")
