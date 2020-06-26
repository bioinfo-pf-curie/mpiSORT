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
    * [Benchmark](#benchmark)
* [Examples](#examples)
    * [Standard](#standard)
    * [Slurm](#slurm)
    * [PBS/Torque](#pbstorque)
* [Performance](#performance)
* [Parallel filesystems](#parallel-filesystems)
* [Algorithm](#algorithm)
* [References](#references)
* [FAQ](#faq)

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

The `-n` options passed to `mpirun` indicates the number of processes to run in parallel (this is basically the number of cores that will be used). For more details on how to choose the number processes, see the [Informatic resources](#informatic-resources) section.

`mpiSORT` requires two mandatory arguments:

* [Input](#input): the path to the SAM file to be sorted
* [Output](#output): the directory in which the results will be written

### Input

A SAM file produced by an aligner (such as [BWA](https://github.com/lh3/bwa)) compliant with the SAM format.

### Options

* `-q INTEGER` filters the reads according to their quality. Reads quality under the threshold are ignored in the sorting results. Default is 0 (all reads are kept) (optional)
* `-p` if the read are paired-end (by default reads are single-end) (optional)
* `-n` sorts the read by their name (but it is not commonly used) (optional)
* `-u` it the input SAM are results of [mpiBWAByChr](https://github.com/bioinfo-pf-curie/mpiBWA) or if there is only one chromosome in the SAM file (optional)

 

### Output

The output consists of gz files:
* one per chromosome (e.g. chr11.gz)
* one for discordant reads (discordant.gz): discordant reads are reads where one pair aligns on a chromosome and the other pair aligns on another chromosome
* one for unmapped reads (unmapped.gz): unmapped reads are reads without coordinates on any chromosome
* when `-u` is set the unmapped and discordant file are prefixed with the chromosome's name

Note that:

1. the discordant.sam file produced by `mpiSORT` is different from the `discordant.sam` produced by `mpiBWAByChr` (see [mpiBWA](https://github.com/bioinfo-pf-curie/mpiBWA) documentation).

2. if you use `mpiSORT` with `-u` option on a chromosome sam file produced by `mpiBWAByChr` (e.g. chr1.sam), `mpiSORT` will produce a file named discordant.sam that will contain the supplementary discordant alignments produced by BWA (i.e the discordant fragments with a flag above 2048). All the primary alignments are sorted. In this case, the `unmapped.sam` will not be produced either by `mpiSORT` as the unmapped reads have already been filtered out by `mpiBWAByChr` in a dedicated file.

3. if you use `mpiSORT` on a sam file produced by `mpiBWA`, the discordant.sam contains all primary and secondary discordant alignments. The unmapped.gz file is generally not empty.


To index the SAM:
`tabix -p sam chr11.gz`

`tabix` is available from [Samtools](http://www.htslib.org/doc/tabix.html)

To uncompress:
`bgzip -d chr11.gz -c > chr11.sam`

`bgzip` is available from [Samtools](http://www.htslib.org/doc/bgzip.html)

## Informatic resources

### Memory

When the SAM file to be sorted contains all the chromosomes, the total memory used during the sorting is around 1.5 times the size of the SAM file.
For example, to sort a 1.3TB SAM file (such as the [NA24631](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/) from [GIAB](https://github.com/genome-in-a-bottle/about_GIAB) which is a 300X Whole Genome (2x150-base) paired reads that we aligned with [mpiBWA](https://github.com/bioinfo-pf-curie/mpiBWA)), 1.7 TB of memory are required and split over 512 MPI workers (i.e. cores) that corresponds to 3.3 GB of memory per core. Note that when the SAM file contains several chromosomes, each chromosome is sorted successively: therefore, the upper memory bound depends on the size of the whole SAM file plus some internal structures plus the size of the biggest chromosome.

NA24631 sample is available here: ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams

To reduce further the memory required you can use as input a SAM file that contains the reads from only one chromosome (as those provided with [mpiBWAByChr](https://github.com/bioinfo-pf-curie/mpiBWA)) and pass it with the option `-u` to `mpiSORT`. The total memory needed in this case is around 2.5 times the size of the individual SAM size. This method also allows a better dispatch of the workload on the computing cluster.

To get a good understanding of the memory management with mpiSORT read with attention the [Benchmark](#benchmark)

### Cpu

Due to the bitonic sorting, the algorithm is optimized for power of 2 number of CPU. Therefore, it is mandatory to set the `-n` parameter of `mpirun` to 2, 4, 8, 16, 32, etc. in order to ensure for optimal performance. For example, `mpirun -n 4 mpiSORT examples/data/HCC1187C_70K_READS.sam ${HOME}/mpiSORTExample`

### Benchmark

This section provides some guidelines to benchmark `mpiSORT` with your infrastructure. It is intended to help the reader to assess what is the best use case and configuration to efficiently benefit from MPI parallelization depending on your computing cluster infrastructure. We strongly recommend that you read carefully this section before running `mpiSORT` on your cluster.

The figures have been obtained from a benchmark using Open MPI 3.1.4 on Intel Skylake architecture.

#### Assess the speed-up with respect to a reference baseline

In this section we are going to compute the baseline walltime for the reference tool `samtools` and assess the speed-up obtained with `mpiSORT`.
We start with 8 cores on a single node and a SAM of 70 GB big containing only 1 chromosome after alignment with [mpiBWAByChr](https://github.com/bioinfo-pf-curie/mpiBWA)).


```
/usr/bin/time -v samtools sort -@ 8 -m 15G -O BAM -o sample70GB.sort.bam -T /tmp/sample70G.prefix.bam sample70GB.sam

Elapsed (wall clock) time (h:mm:ss or m:ss): 42:49.15
Maximum resident set size (kbytes): 68843140

```

The result of the `/usr/bin/time` command shows that `samtools` uses a maximum of ~69 GB of RAM memory (i.e. the resident set size).


```
/usr/bin/time -v mpirun -N 1 -n 8 mpiSort sample70GB.sam -u -q 0

Elapsed (wall clock) time (h:mm:ss or m:ss): 6:46.56
Maximum resident set size (kbytes): 20510832

```

The result of the `/usr/bin/time` command shows that `mpiSORT` uses a maximum of ~ 20 GB of RAM memory per core (i.e. the resident set size).
This means that the total RAM memory needed for mpiSORT is 160 GB = 20 GB * 8 cores.
`mpiSORT` offers a speed-up over 6 with respect to `samtools`.


Now, we increase the number of jobs:


```
/usr/bin/time -v samtools sort -@ 16 -m 10G -O BAM -o sample70GB.sort.bam -T /tmp/sample70G.prefix.bam sample70GB.sam

Elapsed (wall clock) time (h:mm:ss or m:ss): 39:58.98
Maximum resident set size (kbytes): 68859352

```

`samtools` uses the same amout of RAM memory as the previous case.


```
/usr/bin/time -v mpirun -N 1 -n 16 mpiSort sample70GB.sam -u -q 0

Elapsed (wall clock) time (h:mm:ss or m:ss): 3:38.76
Maximum resident set size (kbytes): 10638528

```
Doubling the number of cores decreases by a factor 2 the maximum of RAM memory per core for `mpiSORT` as expected.

Finally, we add a new node to run `mpiSORT`:

```
/usr/bin/time -v mpirun -N 2 -npernode 16 -n 32 mpiSort sample70GB.sam -u -q 0

Elapsed (wall clock) time (h:mm:ss or m:ss): 1:55.31
Maximum resident set size (kbytes): 5003808
```


You normally notice the scalability of `mpiSORT` when the number of cores increases otherwise something is wrong with your setup.

Profiling different sample sizes with different configuration will tell what is the best configuration for dispatching your jobs.
Using this approach can help to decide if you need to use more cores and less GB/core or the contrary. This is very dependent of your computing architecture and your aim.


#### Conclusion

Testing and benchmarking `mpiSORT` with respect to a reference baseline is very important if you want to take advantages of the MPI optimization.
Besides that, there are may MPI parameters you can play with (such as job placement, memory binding and NUMA control).
With this version we recommend you to take a power of 2 MPI jobs in order to avoid memory limitation per jobs.


## Examples

There are many ways to distribute and bind MPI jobs according to your architecture. We provide below several examples to launch MPI jobs in a standard manner or with a job scheduling system such as [Slurm](https://slurm.schedmd.com/sbatch.html) and [PBS/Torque](https://support.adaptivecomputing.com/support/documentation-index/torque-resource-manager-documentation/).

A toy dataset (SAM file) is provided in the [examples/data](../examples/data) folder for testing the program. We test it with from 1 to 8 jobs and 2 nodes with a normal network and file system. You can use this sample to test the program.

### Standard

`mpirun` can be launched in a standard manner without using any job scheduling systems. For example:

`mpirun -n 4 mpiSORT examples/data/HCC1187C_70K_READS.sam ${HOME}/mpiSORTExample -p`

If needed, a file with the server name in `-host` option can be provided to `mpirun`. We invite you to read the `mpirun` documentation for more details.

You can go in the [examples](../examples) directory and execute the [standard.sh](../examples/standard.sh) script to test the program.

### Slurm

In order to submit a job using [Slurm](https://slurm.schedmd.com/sbatch.html), you can write a shell script as follows:

```shell
#! /bin/bash
#SBATCH -J MPISORT_MYSAM_4_JOBS
#SBATCH -N 2                       	# Ask 2 nodes
#SBATCH -n 4                       	# Total number of mpi jobs
#SBATCH -c 1			   	# use 1 core per mpi jobs
#SBATCH --tasks-per-node=2         	# Ask 2 mpi jobs per node
#SBATCH --mem-per-cpu=${MEM}	   	# See Memory ressources
#SBATCH -t 01:00:00
#SBATCH -o STDOUT_FILE.%j.o
#SBATCH -e STDERR_FILE.%j.e

mpirun mpiSORT examples/data/HCC1187C_70K_READS.sam ${HOME}/mpiSORTExample -p -q 0


```

You can go in the [examples](../examples) directory and submit the job with `sbatch` using the [slurm.sh](../examples/slurm.sh) script to test the program.

### PBS/Torque

In order to submit a job using [PBS/Torque](https://support.adaptivecomputing.com/support/documentation-index/torque-resource-manager-documentation/), you can write a shell script as follows:

```shell
#! /bin/bash
#PBS -N MPISORT_MYSAM_4_JOBS
#PBS -l	nodes=2:ppn=2:mem=${MEM}     	# Ask 2 nodes and 2 jobs per node
#PBS -l walltime=24:00:00
#PBS -o STDOUT_FILE.%j.o
#PBS -e STDERR_FILE.%j.e

mpirun mpiSORT examples/data/HCC1187C_70K_READS.sam ${HOME}/mpiSORTExample -p -q 0

```

You can go in the [examples](../examples) directory and submit the job with `qsub` command using the [pbs.sh](../examples/pbs.sh) script to test the program.

## Parallel filesystems

As the `mpiSORT` program uses MPI functions for reading and writing you can take advantage of a parallel file system to tackle the IOs bottleneck and speed-up the sorting of a SAM file. Using a parallel filesystem such as [Lustre](http://lustre.org/) or [BeeGFS](https://www.beegfs.io/) will be mandatory as long as the SAM file to sort is bigger and bigger.


Note that for advanced users, it is even possible to fine tune the source code such that `mpiSORT` can optimally used the settings of the parallel filesystem on which you will run the program. Refer to the [Advanced](ADVANCED.md) guidelines to do so.


We don't recommend to use MPI with [NFS](https://en.wikipedia.org/wiki/Network_File_System) (it works but it does not scale very well).

## Performance

Obviously, the performance of `mpiSORT` depends on the computing infrastructure. Using the computing cluster provided by the [TGCC France Génomique](https://www.france-genomique.org/plateformes-et-equipements/plateforme-tgcc-arpajon/) (Broadwell architecture), we obtained the following performance sorting the [NA24631](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/) from [GIAB](https://github.com/genome-in-a-bottle/about_GIAB) which is a 300X Whole Genome:


* **20 minutes** with 512 jobs and with an efficiency of 70% (30% time is spent in IO) with a [Lustre](http://lustre.org/) configuration :
    * `lfs setstripe -c 12 -S 4m` (for input)
    * `lfs setstripe -c 12 -S 256m` (for output)

* **10 minutes** with 1024 jobs with an efficiency of 60% (40% time is spent in IO) with a Lustre configuration :
    * `lfs setstripe -c 12 -S 256m` (for input)
    * `lfs setstripe -c 12 -S 256m` (for output)

Because of the development design the program is optimized for HPC architecture. This program runs better on low latency network and parallel file system. 

## Algorithm

Sorting a file is all about IO's and shuffling (or movements) of data. Therefore, we developed a program that capitalizes on parallel and distributed file system in order to overcome the bottlenecks encountered with traditional tools to sort large sequencing data. Our approach relies on message passing interface paradigm (MPI) and distributed memory available on high computing performance architecture.

The program we implemented in based on the major components: the bitonic-sort, the shuffling of the data and the distributed cache.

The [bitonic sort](https://en.wikipedia.org/wiki/Bitonic_sorter) is a real parallel sorting algorithm that works on parallel architectures. The complexity of the bitonic is of n(log(n))^2. The bitonic sorter has been developed using MPI message passing primitives and is inspired from the book of [Peter S. Pacheco "Parallel programming with MPI".](https://www.cs.usfca.edu/~peter/ppmpi/).


![Bitonic sort](https://upload.wikimedia.org/wikipedia/commons/b/bd/BitonicSort1.svg)
[bitonic sort](https://en.wikipedia.org/wiki/Bitonic_sorter)

The shuffling of the data is done through the Bruck method. This method has the advantage of avoiding the shuffle bottleneck (The All2all). Bruck is a log(N) method and scale very well for distributed architectures.

For more details, see the [References](#references) section.

## References

Original Bruck algorithm:

* Bruck et al. (1997) [Efficient Algorithms for All-to-All Communications in Multiport Message-Passing Systems](https://dl.acm.org/doi/10.1109/71.642949). EEE Transactions on Parallel and Distributed Systems, 1997 ([pdf](http://authors.library.caltech.edu/12348/1/BRUieeetpds97.pdf)).

Modified versions of the Bruck algorithm:

* Jesper Larsson Träff, Antoine  Rougier and Sascha  Hunold, [Implementing a Classic:Zero-copy All-to-all Communication with MPI Datatypes](https://dl.acm.org/doi/10.1145/2597652.2597662). ICS '14: Proceedings of the 28th ACM international conference on Supercomputing, 2014 ([pdf](http://hunoldscience.net/paper/classical_sahu_2014.pdf)).

* Rajeev Thakur, Rolf  Rabenseifner and William Gropp, [Optimization of Collective Communication Operations in MPICH](https://dl.acm.org/doi/10.1177/1094342005051521). International Journal of High Performance Computing Applications, 2005 ([pdf](http://www.mcs.anl.gov/~thakur/papers/ijhpca-coll.pdf)).

Bitonic sort:

* Ananth Grama, Anshul Gupta, George Karypis, Vipin Kumar, [Introduction to Parallel Computing](https://www-users.cs.umn.edu/~karypis/parbook/). Addison-Wesley Edition, 2003 ([pdf](http://srmcse.weebly.com/uploads/8/9/0/9/8909020/introduction_to_parallel_computing_second_edition-ananth_grama..pdf)).

* Kenneth E. Batcher, Sorting networks and their applications. American Federation of Information Processing Societies: AFIPS Conference Proceedings: 1968 Spring Joint Computer Conference, Atlantic City, NJ, USA ([pdf](Sorting networks and their applications))

Description of the SAM format:

* Li H. et al, [The sequence alignment/map format and SAMtools](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/). Bioinformatics, 2009.


* [samtools](https://github.com/lh3/samtools) source code. The `bgzf.c` program was used and included in the `mpiSORT` source code. [samtools](https://github.com/lh3/samtools) is provided under *The MIT License,  Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology*.

Presentations about our program:

* Journées nationales du DEVeloppement logiciel, 2015 [pdf](http://devlog.cnrs.fr/_media/jdev2015/poster_jdev2015_institut_curie_hpc_sequencage_jarlier.pdf?id=jdev2015%3Aposters&cache=cache)
* OpenSFS conference Lustre LAD, 2016 [pdf](http://www.eofs.eu/_media/events/lad16/03_speedup_whole_genome_analysis_jarlier.pdf)
* Journées nationales du DEVeloppement logiciel, 2017 [pdf](http://devlog.cnrs.fr/_media/jdev2017/poster_jdev2017_hpcngs_frederic_jarlier.pdf?id=jdev2017%3Aposters&cache=cache)

## FAQ

### Is it a mandatory to have a power of 2 keys (reads) in my SAM files?

No it is not a mandatory. It is best case but in practice and most of the time it is impossible (you can think of the FFTW tool). The sorting ask for a power of 2 keys to work so internally we fill the coordinates vector to fit that power of 2. You may think this will occur an important overhead, but in fact not to much, around 9% in time with N=2^20 and N=2^21  (according to your computing infrastructure). Why? Because the algorithm works in O(log²) time complexity. So whether the number of keys is between 2^N and 2^(N+1) the lower and upper bounds sorting time are set and your sorting time is predictable (or real time). In conclusion each time the keys increase with a factor of 2 the time increases by (N/(N+1))^2. The counterpart is you will need more resources.  
   
### Is it a mandatory to have a power of 2 number of CPU to use mpiSORT?

With this actual version yes. If you really want to play with non power of 2 CPU no problem, in mpiSort.c comment from the lines 220 to 230 and recompile the source. Using non power of 2 CPU reduce the number of CPU resources needed but induce a memory overhead for the rank 0 job (some MB). Try it and if it works then use it but be aware of the memory bound (see [Where does this memory bounds comes from?](#where-does-this-memory-bounds-comes-from))

### Where does this memory bounds comes from?

The memory bounds have been released for power of 2 cores but is still in place if you decided to use non-power of 2 cores. The work around is to send buffer by packets of 1GB during the shuffle of the reads.   

### Is the sorting stable?

No the bitonic sort is not stable. It means the order of the reads with same coordinates may vary according to parallelization. Does it really matter? well its a debate because the order of the reads in the fastq files is already a random process once you use parallelization during alignment. This is also why you cannot use md5 to test reproducibility with sorted files unless you take the same number of parallel process. Stabilization could be done with adding a second key during the sorting, this will be an option in the future release.     

