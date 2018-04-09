#!/usr/bin/env python

import argparse
import subprocess as sp
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
    if int(genomeSizeMb) != 0:
        genomeSize = int(genomeSizeMb)*1000000
    else:
        genomeSize = 1000000
    desiredAmountOfData = genomeSize*int(desiredCoverage)
    sequenceLengths = np.array(sequenceLengths)
    sequenceLengths = np.sort(sequenceLengths)[::-1]
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
    if idx > len(sequenceLengths):
        print "You do not have enough sequence data to attain " + str(desiredCoverage) + "X coverage. "
    elif idx <= len(sequenceLengths):
        cutoff = sequenceLengths[idx]
        print "Retain sequences larger than "+str(cutoff)+" to achieve "+str(desiredCoverage)+"X coverage. "

    # plot accumulation curve
    
    fig = plt.figure()
    plt.scatter(fewerX,fewerY, c="#d62728")
    plt.xlim(max(x)*-0.05, max(x)+max(x)*0.05)
    plt.ylim(max(y)*-0.1, max(y)+max(y)*0.1)
    yTickLabels = []
    if max(y) >= genomeSize*10:
        yArrayOfTicks = np.arange(0, max(y), genomeSize*10)
    if max(y) < genomeSize*10:
        yArrayOfTicks = np.arange(0, max(y), genomeSize*1)
    for tick in np.nditer(yArrayOfTicks):
        if int(genomeSizeMb) == 0:
            oneLabel = str(tick/(genomeSize*0.1))
            plt.ylabel("Data (Mb)")
        else:
    	    oneLabel = str(tick/genomeSize) + "x"
            plt.ylabel("Genome coverage")
	yTickLabels.append(oneLabel)
    plt.yticks(yArrayOfTicks, yTickLabels, fontsize=10)
    plt.title("Accumulation curve")
    # plt.show()
    if int(desiredCoverage) == 1 or idx > len(sequenceLengths):
        plt.xlabel("Nth read")
        xArrayOfTicks = np.arange(0, max(x), len(x)*0.1)
        plt.xticks(xArrayOfTicks, fontsize=10, rotation='vertical')
        plt.savefig(outputPrefix+"_accumulationCurve.pdf", format='pdf')
    elif int(desiredCoverage) > 1:
        plt.plot([max(x)*-0.05, idx], [desiredAmountOfData, desiredAmountOfData], color='r', linestyle='-')
        plt.plot([idx, idx], [max(y)*-0.1, desiredAmountOfData], color='r', linestyle='-')
        plt.xticks([])
        cutoffText = str(cutoff) + "bp"
        plt.text(idx, max(y)*-0.1275, cutoffText, rotation='vertical', fontsize=8, ha='center', ma='right')
        plt.savefig(outputPrefix+"_accumulationCurve_"+str(desiredCoverage)+"X.pdf", format='pdf')

    plt.close(fig)

    # plot read length histogram

    fig = plt.figure()
    plt.bar(range(len(myBins)), countsInBins, color="#d62728", align='edge')
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

parser = argparse.ArgumentParser(description="Visualize long read data from sequence files (.fasta or .fastq, zipped or unzipped)")
parser.add_argument("sequenceFile", type=argparse.FileType('r'), help=".fastq or .fastq.gz file containing all your sequences.")
parser.add_argument("-out", "--outputPrefix", help="Prefix for all output files.", type=str, default="out")
parser.add_argument("-gs", "--genomeSize", help="The genome size in Mb", type=int, default=0)
parser.add_argument("-cov", "--coverage", help="Your desired genome coverage", type=int, default=1)
args = parser.parse_args()

if args.sequenceFile.name.endswith(".gz"):
    command = "gunzip -c " + args.sequenceFile.name + " >" + args.sequenceFile.name.strip(".gz")
    sp.check_output(command, shell=True)
    tempfastq = args.sequenceFile.name.strip(".gz")
    nanoporeplots(tempfastq, args.outputPrefix, args.genomeSize, args.coverage)
    sp.check_output("rm " + tempfastq, shell=True)
elif args.sequenceFile.name.endswith("q"):
    nanoporeplots(args.sequenceFile.name, args.outputPrefix, args.genomeSize, args.coverage)

