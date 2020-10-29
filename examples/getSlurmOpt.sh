#! /usr/bin/python

import argparse
import math

parser = argparse.ArgumentParser(description = 'Estimate informatic resources for mpiSORT')

parser.add_argument('-c','--cores', dest = 'cores', help = 'Number of cores in a node', type = int, required = True)
parser.add_argument('-m','--memory', dest = 'memory', help = 'Maximum RAM memory available in a node in GB', type = float, required = True)
parser.add_argument('-s','--size', dest = 'size', help = 'Size of the SAM file in GB', type = float, required = True)
parser.add_argument('-t','--type', dest = 'type', help = 'Type of SAM (uniq or all)', type = str, required = True)

args = parser.parse_args()

### When the FASTQ file contains all the chromosomes,
### the memory needed is the size of the FASTQ file multiplied by scaleFactorAllChr
scaleFactorAllChr = 1.5

### When the FASTQ file contains only one chromosome,
### the memory needed is the size of the FASTQ file multiplied by scaleFactorOneChr
scaleFactorOneChr = 3.5


### Functions
def computeInformaticResources(cores, size, memory, scaleFactor):
    
    requiredMemory = size * scaleFactor
    requiredNodes = math.ceil(requiredMemory / memory)
    memoryByCore = memory / cores
    requiredCores = math.ceil(requiredMemory / memoryByCore)
    power2Cores = math.floor(math.log(math.ceil(requiredCores)) / math.log(2))
    power2Cores = 2 ** power2Cores
    power2MemoryPerCore = round(requiredMemory / power2Cores, 1)
    power2CoresByNode = math.ceil(power2Cores / requiredNodes)

    return requiredCores, requiredMemory, requiredNodes, power2Cores, power2MemoryPerCore, power2CoresByNode


def printResources(resources, message):
    print("\n")
    print(message + ", the informatic resources required are:")
    print("\tMemory: {} GB".format(round(resources[1], 1)))
    print("\tNumber of cores: {}".format(resources[0]))
    print("\tNumber of nodes: {}".format(resources[2]))
    print("\tNumber of cores (with power of 2 constraint): {}".format(resources[3]))
    print("\tMemory per of core : {}".format(resources[4]))
    print("\tNumber of cores by node (with power of 2 constraint): {}".format(resources[5]))

def printCommandLine(resources):
    print("\t-N {} -n {} -c 1 --tasks-per-node={} --mem-per-cpu={}GB\t".format(int(resources[2]), int(resources[3]), int(resources[5]),int(resources[4]) + 1))


### Compute informatic resources
if type in ['all']:
	resourceAllChr = computeInformaticResources(args.cores, args.size, args.memory, scaleFactorAllChr)
        printCommandLine(resourceAllChr)
else:
	resourceOneChr = computeInformaticResources(args.cores, args.size, args.memory, scaleFactorOneChr)
        printCommandLine(resourceOneChr)





 
