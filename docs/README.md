# Documentation

* [Installation](#installation)
* [Prerequisites](#prerequisites)
* [Usage](#usage)
    * [Input](#input)
    * [Options](#options)
    * [Output](#output)
* [Informatic resources](#informatic-resources)
    * [Memory](#memory)
    * [Cpu](#cpu)
* [Examples](#examples)
    * [Standard](#standard)
    * [Slurm](#slurm)
    * [PBS/Torque](#pbstorque)
* [Performance](#performance)
* [Filesystems](#filesystems)
* [Algorithm](#algorithm)
* [References](#references)


## Installation

Follow the [Installation](INSTALL.md) guidelines to compile and install `mpiSORT`.

## Prerequisites

As `mpiSORT` relies on the Message Passing Interface (MPI) standard, `mpirun` must be available to run the program. Several MPI implementations exist
such as [mpich](https://www.mpich.org/), [open-mpi](https://www.open-mpi.org/) or [Intel® MPI Library](https://software.intel.com/en-us/mpi-library). The `mpirun` program must be available in your PATH.

### Install mpirun on CentOS

* for [open-mpi](https://www.open-mpi.org/): `sudo yum install openmpi openmpi-devel`
* for [mpich](https://www.mpich.org/): `sudo yum install mpich`

### Install mpirun on Ubuntu

* for [open-mpi](https://www.open-mpi.org/): `sudo apt-get install openmpi-bin`
* for [mpich](https://www.mpich.org/): `sudo apt-get install mpich`



## Usage

`mpiSORT` is executed with the `mpirun` program, for example:

`mpirun -n 4 mpiSORT examples/data/HCC1187C_70K_READS.sam ${HOME}/mpiSORTExample`

The `-n` options passed to `mpirun` indicates the number of processes to run in parallel (this is basically the number of cores that will be used). For for details on how to choose the number processes, see the [Informatic resources](#informatic-resources) section.

`mpiSORT` requires two mandatory arguments:

* [Input](#input): the path to the SAM file to be sorted
* [Output](#output): the directory in which the results will be written

### Input

A SAM file produced by an aligner (such as [BWA](https://github.com/lh3/bwa)) with paired reads and compliant with the SAM format. Only paired paired.

### Options

* `-q INTEGER` filters the reads according to their quality. Reads quality under the threshold are ignored in the sorting results. Default is 0 (all reads are kept).
* `-n` sorts the read by their name (but it is not commonly used).

### Output

The output consists of gz files:
* one per chromosome (e.g. chr11.gz)
* one for discordant reads (discordant.gz)
* one for unmapped reads (unmapped.gz)

To index the SAM:
`tabix -p sam chr11.gz`

`tabix` is available from [Samtools](http://www.htslib.org/doc/tabix.html)

To uncompress:
`bgzip -d chr11.gz -c > chr11.sam`

`bgzip` is available from [Samtools](http://www.htslib.org/doc/bgzip.html)

## Informatic resources

### Memory

The total memory used during the sorting is around one and a half the size of the SAM file.
For example, to sort a 1.3TB SAM file (such as the [NA24631](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/) from [GIAB](https://github.com/genome-in-a-bottle/about_GIAB) which is a 300X Whole Genome (2x150-base) paired reads that we aligned with [mpiBWA](https://github.com/bioinfo-pf-curie/mpiBWA)), 1.7 TB of memory are required and splitted over 512 MPI workers (i.e. cores) that corresponds to makes 3.3 Gb of memory per core.

NA24631 sample is available here: ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams

### Cpu

Due to the bitonic sorting, the algorithm is optimized for power of 2 number of CPU. Therefore, it is recommended to set the `-n` parameter of `mpirun` to 2, 4, 8, 16, 32, etc. in order to ensure for optimal performance. For example, `mpirun -n 4 mpiSORT examples/data/HCC1187C_70K_READS.sam ${HOME}/mpiSORTExample`


However, the  `-n` parameter can be set to any other value but extra  MPI communications will be added to fit power of 2 required by the bitonic algorithm. In this case, additonal memory is needed for the rank 0 worker. This rank is responsible for collecting and dispatching the data before and after bitonic.


## Examples

There are many ways to distribute and bind MPI jobs according to your architecture. We provide below several examples to launch MPI jobs in a standard manner or with a job scheduling system such as [Slurm](https://slurm.schedmd.com/sbatch.html) and [PBS/Torque](https://support.adaptivecomputing.com/support/documentation-index/torque-resource-manager-documentation/).

A toy dataset (SAM file) is provided in the [examples/data](../examples/data) folder for testing the program. We test it with from 1 to 8 jobs and 2 nodes with a normal network and file system. You can use this sample to test the program.

### Standard

`mpirun` can be launched in a standard manner without using any job scheduling systems. For example:

`mpirun -n 4 mpiSORT examples/data/HCC1187C_70K_READS.sam ${HOME}/mpiSORTExample`

If needed, a file with the server name in `-host` option can be provided to `mpirun`. We invite you to read the `mpirun` documentation for more details.


### Slurm

```shell
#! /bin/bash
#SBATCH -J MPISORT_MYSAM_64_CORES
#SBATCH -N 2                       # Ask 2 nodes
#SBATCH -n 4                      # Total number of cores
#SBATCH -c 1
#SBATCH --tasks-per-node=2         # Ask 2 cores per node
#SBATCH -t 01:00:00
#SBATCH -o STDOUT_FILE.%j.o
#SBATCH -e STDERR_FILE.%j.e

mpirun mpiSORT examples/data/HCC1187C_70K_READS.sam ${HOME}/mpiSORTExample -q 0

```

### PBS/Torque


## Filesystems


As the programs use MPI fonctions for reading and writing you can take advantage of a parallel file system. To speed-up reading and writing you can set the striping of your data. The striping tells the number of file servers you want to use and the size of data blocks on each server. You can compare the striping of the file with the mapping process used in Hadoop. This is the way your data are distributed among the servers of your file system. This kind of optimization accelerate drastically the IO operations.

Standard software (Samtools, Sambamba, Picard,... ) don't take into account the underlying distributed file system and low latency interconnexion when MPI does.

### Parallel file system configuration:

MPI improves IOs by means of parallelization.

On parallel file system like Lustre, GPFS one way of accelerating IO is to stripe them across servers.
You can choose the number of servers with the striping factor and the size of chunks with the striping size. 

Before running and compile you must tell the programm how the input data is striped on Lustre and how to stripe the output data.
The parallel writing and reading are done via the MPI finfo structure in the first line of mpiSort.c.

For reading and set the Lustre buffer size.
To do that edit the code of mpiSort.c and change the parameters in the header part.
After tuning parameters recompile the application.  

If you are familiar with MPI IO operations you can also test different commands collective, double buffered, data sieving.
In file write.c in the function writeSam, writeSam_unmapped, writeSam_discordant. 

The default parameters of the programm are unharmed for other filesystem. 

At TGCC with a striping factor of 128 (number of OSSs servers) and the a striping size of 2.5 GB, a SAM file of 1.2TB could be loaded in memory in less than 1 minute.  

### Lustre

The section of the code you can modify in MPI info according to your Lustre configuration:

For reading part (mpiSort.c): 
line 98 => 102 and 261 => 280  

For writing part (in write.c):  
line 1406 => 1412  
line 3516 => 3521  

The parameters are harmless for other filesystem.  

!!!We recommand to test different parameters before setting them once for all.  

### Network File system configuration:

We don't recommend to use MPI with NFS (it works but it does not scale very well). 

If you run the programm on NFS take a power of 2 number of jobs. 

For NFS users there is a bug in openMPI. 
https://www.open-mpi.org/community/lists/users/2016/06/29434.php

a patch is available here:
https://trac.mpich.org/projects/mpich/attachment/ticket/2338/ADIOI_NFS_ReadStrided.patch

To improve TCP communications over openMPI : 
https://www.open-mpi.org/faq/?category=tcp

3) Algorithm 
12) Parallel file system configuration 
13) Lustre optimization
14) Network file system configuration 

## Performance

Obviously, the performance of `mpiSORT` depends on the computing infrastruture. Using the computing cluster provided by the [TGCC France Génomique](https://www.france-genomique.org/plateformes-et-equipements/plateforme-tgcc-arpajon/) (Broadwell architecture), we obtained the following performance sorting the [NA24631](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/) from [GIAB](https://github.com/genome-in-a-bottle/about_GIAB) which is a 300X Whole Genome:


* **20 minutes** with 512 jobs and with an efficiency of 70% (30% time is spent in IO) with a [Lustre](http://lustre.org/) configuration :
    * `lfs setstripe -c 12 -S 4m` (for input)
    * `lfs setstripe -c 12 -S 256m` (for output)

* **10 minutes** with 1024 jobs with an efficiency of 60% (40% time is spent in IO) with a Lustre configuration :
    * `lfs setstripe -c 12 -S 256m` (for input)
    * `lfs setstripe -c 12 -S 256m` (for output)

Because of the development design the programm is optimized for HPC architecture. This programm runs better on low latency network and parallel file system. 

## Algorithm

Sorting a file is all about IO's and shuffling (or movements) of data. Therefore, we developed a program that capitalizes on parallel and distributed file system in order to overcome the bottlenecks encountered with traditionnal tools to sort large sequencing data. Our approach relies on message passing interface paradigm (MPI) and distributed memory available on high computing performance architecture.

The program we implemented in based on the major components: the bitonic-sort, the shuffling of the data and the distributed cache.

The [bitonic sort](https://en.wikipedia.org/wiki/Bitonic_sorter) is a real parallel sorting algorithm that works on parallel architectures. The complexity of the bitonic is of (log(n))^2 instead of nlog(n) with the parallel merge-sort. The bitonic sorter has been developped using MPI message passing primitives and is inspired from the book of [Peter S. Pacheco "Parallel programming with MPI".](https://www.cs.usfca.edu/~peter/ppmpi/)

The shuffling of the data is done through the Bruck method. This method has the advantage of avoiding the shuffle bottleneck (The All2all). Bruck is a log(N) method and scale very well for distributed architectures.

For more details, see the [References](#references) section.

## References

Original Bruck algorithm:

* Bruck et al. (1997) [Efficient Algorithms for All-to-All Communications in Multiport Message-Passing Systems](https://dl.acm.org/doi/10.1109/71.642949). EEE Transactions on Parallel and Distributed Systems, 1997 ([pdf](http://authors.library.caltech.edu/12348/1/BRUieeetpds97.pdf)).

Modified versions of the Bruch algorithm:

* Jesper Larsson Träff, Antoine  Rougier and Sascha  Hunold, [Implementing a Classic:Zero-copy All-to-all Communication with MPI Datatypes](https://dl.acm.org/doi/10.1145/2597652.2597662). ICS '14: Proceedings of the 28th ACM international conference on Supercomputing, 2014 ([pdf](http://hunoldscience.net/paper/classical_sahu_2014.pdf)).

* Rajeev Thakur, Rolf  Rabenseifner and William Gropp, [Optimization of Collective Communication Operations in MPICH](https://dl.acm.org/doi/10.1177/1094342005051521). International Journal of High Performance Computing Applications, 2005 ([pdf](http://www.mcs.anl.gov/~thakur/papers/ijhpca-coll.pdf)).

Description of the SAM format:

* Li H. et al, [The sequence alignment/map format and SAMtools](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/). Bioinformatics, 2009.


* [samtools](https://github.com/lh3/samtools) source code. The `bgzf.c` program was used and included in the `mpiSORT` source code. [samtools](https://github.com/lh3/samtools) is provided under *The MIT License,  Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology*.

Presentation about our program:

* Journées nationales du DEVeloppement logiciel, 2015 [pdf](http://devlog.cnrs.fr/_media/jdev2015/poster_jdev2015_institut_curie_hpc_sequencage_jarlier.pdf?id=jdev2015%3Aposters&cache=cache)
* OpenSFS conference Lustre LAD, 2016 [pdf](http://www.eofs.eu/_media/events/lad16/03_speedup_whole_genome_analysis_jarlier.pdf)
* Journées nationales du DEVeloppement logiciel, 2017 [pdf](http://devlog.cnrs.fr/_media/jdev2017/poster_jdev2017_hpcngs_frederic_jarlier.pdf?id=jdev2017%3Aposters&cache=cache)
