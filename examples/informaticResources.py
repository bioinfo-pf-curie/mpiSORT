#! /usr/bin/python

import argparse
import math

parser = argparse.ArgumentParser(description = 'Estimate informatic resources for mpiSORT')

parser.add_argument('-c','--cores', dest = 'cores', help = 'Number of cores in a node', type = int, required = True)
parser.add_argument('-m','--memory', dest = 'memory', help = 'Maximum RAM memory available in a node in GB', type = float, required = True)
parser.add_argument('-s','--size', dest = 'size', help = 'Size of the SAM file in GB', type = float, required = True)

args = parser.parse_args()

### When the FASTQ file contains all the chromosomes,
### the memory needed is the size of the FASTQ file multiplied by scaleFactorAllChr
scaleFactorAllChr = 1.5

### When the FASTQ file contains only one chromosome,
### the memory needed is the size of the FASTQ file multiplied by scaleFactorOneChr
scaleFactorOneChr = 2.5


### Functions
def computeInformaticResources(cores, size, memory, scaleFactor):
    
    requiredMemory = size * scaleFactor
    requiredNodes = requiredMemory / memory
    memoryByCore = memory / cores
    requiredCores = requiredMemory / memoryByCore

    return requiredCores, requiredMemory, requiredNodes


def printResources(resources, message):
    print("\n")
    print(message + ", the informatic resources required are:")
    print("\tMemory: {} GB".format(round(resources[1], 1)))
    print("\tNumber of nodes: {}".format(math.ceil(resources[2])))
    print("\tNumber of cores: {}".format(math.ceil(resources[0])))


### Compute informatic resources
resourceAllChr = computeInformaticResources(args.cores, args.size, args.memory, scaleFactorAllChr)
resourceOneChr = computeInformaticResources(args.cores, args.size, args.memory, scaleFactorOneChr)

### Print the setting
print("\n")
print("Your setting is:")
print("\tA node has {} cores".format(args.cores))
print("\tA node has {} GB of RAM memory".format(args.memory))
print("\tThe size of the SAM file is {} GB".format(args.size))

### Print informatic resources
printResources(resourceAllChr, "If your FASTQ file contains all the chromosomes")
printResources(resourceOneChr, "If your FASTQ file contains only one chromosome")
