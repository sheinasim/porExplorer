#!/usr/bin/env python

import subprocess as sp
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import math

# usage
def usage():
    print("\n-usage: "+sys.argv[0]+" <fastq> <outputPrefix> <genome size in Mb (write 100 Mb genome as 100)> <desired coverage>")  
    sys.exit()

def calculateSumsInBins(arrayOfBins, arrayOfSequenceLengths):
    listOfSums = []
    for myBin in np.nditer(arrayOfBins):
        binSum = np.sum(arrayOfSequenceLengths[np.where((arrayOfSequenceLengths < myBin) & (arrayOfSequenceLengths > (myBin-1000)))])
        listOfSums.append(binSum)
    return listOfSums

def countSeqsInBins(arrayOfBins, arrayOfSequenceLengths):
    listOfCounts = []
    for myBin in np.nditer(arrayOfBins):
        binCounts = len(arrayOfSequenceLengths[np.where((arrayOfSequenceLengths < myBin) & (arrayOfSequenceLengths > (myBin-1000)))])
        listOfCounts.append(binCounts)
    return listOfCounts

def nanoporeplots(inFastq, outputPrefix="output", genomeSizeMb=0, desiredCoverage=1):
    sequenceLengths = []
    with open(inFastq, 'r') as infile:
    	for line in infile:
	    line = line.strip()
	    if line.startswith(("A", "C", "G", "T", "N")):
	        sequenceLengths.append(len(line))
	    else:
	        pass
    genomeSize = int(genomeSizeMb)*1000000
    desiredAmountOfData = genomeSize*int(desiredCoverage)
    sequenceLengths = np.array(sequenceLengths)
    sequenceLengths = np.sort(sequenceLengths)[::-1]
    # sequenceLengths = np.flipud(sequenceLengths)
    cumulativeSum = np.cumsum(sequenceLengths)
    x = np.arange(len(cumulativeSum))
    y = cumulativeSum
    fewerX = x[0:len(x):len(x)*0.01]
    fewerY = y[0:len(y):len(y)*0.01]
    myBins = np.arange(1000, 70000, 1000, dtype=float) 
    sumsInBins = np.array(calculateSumsInBins(myBins, sequenceLengths), dtype=float)
    proportionsInBins = sumsInBins/sum(sequenceLengths)
    countsInBins = np.array(countSeqsInBins(myBins, sequenceLengths), dtype=float)
    minimumSum = np.array(cumulativeSum[cumulativeSum < desiredAmountOfData])
    idx = len(minimumSum)+1
    cutoff = sequenceLengths[idx]

    # plot accumulation curve
    
    print "Retain sequences larger than "+str(cutoff)+" to achieve "+str(desiredCoverage)+"X coverage. "
    fig = plt.figure()
    plt.scatter(fewerX,fewerY, c="#328AFF")
    plt.xlim(max(x)*-0.05, max(x)+max(x)*0.05)
    plt.ylim(max(y)*-0.1, max(y)+max(y)*0.1)
    plt.plot([max(x)*-0.05, idx], [desiredAmountOfData, desiredAmountOfData], color='r', linestyle='-')
    plt.plot([idx, idx], [max(y)*-0.1, desiredAmountOfData], color='r', linestyle='-')
    yTickLabels = []
    yArrayOfTicks = np.arange(0, max(y), genomeSize*10)
    for tick in np.nditer(yArrayOfTicks):
    	oneLabel = str(tick/genomeSize) + "x"
	yTickLabels.append(oneLabel)
    plt.yticks(yArrayOfTicks, yTickLabels)
    plt.xticks([])
    plt.ylabel("Genome coverage")
    cutoffText = str(cutoff) + "bp"
    plt.text(idx, max(y)*-0.1275, cutoffText, rotation='vertical', fontsize=8, ha='center', ma='right')
    plt.title("Accumulation curve")
    # plt.show()
    plt.savefig(outputPrefix+"_accumulationCurve_"+str(desiredCoverage)+"X.pdf", format='pdf')
    plt.close(fig)

    # plot read length histogram

    fig = plt.figure()
    plt.bar(range(len(myBins)), countsInBins, color="#328AFF", align='edge')
    plt.title("Histogram of sequence lengths")
    plt.ylabel("Number of reads")
    plt.xlabel("Length (Kb)")
    # plt.show()
    plt.savefig(outputPrefix+"_rawHist.pdf", format='pdf')
    plt.close(fig)
    
    # plot read length histogram as proportion of dataset for real
    fig = plt.figure()
    plt.bar(range(len(myBins)), proportionsInBins, color="#328AFF", align='edge')
    plt.title("Histogram of sequence lengths as a proportion of data")
    plt.ylabel("Proportion of reads")
    plt.xlabel("Length (Kb)")
    # plt.show()
    plt.savefig(outputPrefix+"_proportionHist.pdf", format='pdf')
    plt.close(fig)

if len(sys.argv) < 2:
    usage()
else:
    if sys.argv[1].endswith(".gz"):
        command = 'gunzip -c ' + sys.argv[1] + ' >' + sys.argv[1].strip('.gz')
        sp.check_output(command, shell=True)
        tempfastq = sys.argv[1].strip(".gz")
        nanoporeplots(tempfastq, sys.argv[2], sys.argv[3], sys.argv[4])
        #print tempfastq
        sp.check_output('rm ' + tempfastq, shell=True)
    else:
        nanoporeplots(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

