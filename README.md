Objective
---------

Sorting big NGS data file Version 1.0. 

Sections:
-------

[release notes] 1) Release notes <br />
[installation] 2) Installation <br />
3)  Algorithm <br />
4) Requirements <br />
5) Memory usage <br />
6) Inputs <br />
8) Outputs <br />
9) Compiler <br />
10) Sample test
11) Configuration <br />
11) Lustre configuration <br />
12) File system management <br />
13) Results <br />
14) Job rank <br />
15) Example of code <br />
16) Options <br />
17) Improvements and future work <br />

==============================================================================


# release notes  
-------------

Release 1.0 from 06/10/2016

1) The previous version didn't sort the offset destination before the shuffle. This bug is fixed. <br />
2) New packaging with autools for distribution. <br />

Release 0.9 from 29/07/2016

1) Bug fix when data sample is too small and some jobs have nothing to sort <br />
2) Bug fix when jobs under 6 <br />
3) We add new option -n for sorting by name <br />
4) the parallel merge sort has been replaced with bitonic sorting (25% percent gain in speed-up) <br />
5) Refactoring and cleaning of the code <br />

## installation
-----------

You need automake 1.15 for the installation. <br />
You can install automake and autoconf in differents directories and export the path like this: <br />
export PATH=../automake-1.15/bin:../autoconf-2.69/bin:$PATH <br />

Download from git. In the folder mpiSORT type: <br />
./configure && make install && make <br />

or for distribution: <br />
make dist  <br />
tar xzf .tar.gz  <br />
cd mpisort-1.0  <br />
./configure && make install && make<br />

for passing mpi path: <br />
./configure CC=mpi_bin_path  <br />
add --prefix in configure if you need  <br />

3) Algorithm
----------

Sorting a file is all about IO's and shuffling of data. 

We have developed a real parallel and distributed file system aware program to overcome some issues encounter with traditionnal tools like Samtools, Sambamba, Picard. We propose a novel approach based on distributed memory computer. <br />

There are several aspects in this sorting algorithm that are important: the bitonic-sort, the shuffling of the data and the distributed cache buffering. <br />

The parallel merge-sort has been replaced with a bitonic merge-sort. The bitonic sort is a real parallel sorting algorithm (https://en.wikipedia.org/wiki/Bitonic_sorter). The complexity of the bitonic is of (log(n))^2 instead of nlog(n) with the parallel merge-sort. The bitonic sorter has been developped using MPI message passing primitive. This is why the program runs faster with low latency network.<br />

The shuffing of the data is done through the Bruck method. This method has the advantage of avoiding the shuffle bottleneck (The All2all). Bruck is a log(N) method and scale very well for distributed architectures. <br /> 

The program makes an intensive use of IO reads cache buffering at the client level. The cache size depends of the stripping of your data on the distributed file system. The striping tells the number of file servers you want to use and the size of data blocks on each server. You can compare the striping of the file with the mapping process use in Hadoop. This is the way your data are distributed among the servers of your file system. This kind of optimizations accelerate drastically the IO operations and file management.<br />

Ordinary softwares (Samtools, Sambamba, Picard,... ) doesn't take into account the underlying distributed file system and low latency interconnexion when MPI does. <br />

4) Requirements:
-------------

For small data samples a normal network and a network file system can do the job. <br />

The presented version of the program has been tested on France Genomic cluster of TGCC (Très Grand Centre de Calcul) of CEA (Bruyeres le Chatel, France). <br />

Because of the development design and feature program have been optimized for HPC architecture.

!!!!! THIS PROGRAM RUN BETTER ON A LOW LATENCY NETWORK !!!! <br />
!!!!! THIS PROGRAM RUN BETTER ON A PARALLEL FILE SYSTEM !!!!

This programm is intended to run on supercomputer architectures equiped with Infiniband network and parallel 
filesystem such as Lustre, GPFS,... 

Contact us if you need information.

5) Memory
---------

For a 1 TB SAM file (paired 100pb, 100X coverage) less than 4gb per jobs is needed over 512 cores. 
A peak of 10gb can be observed on jobs master rank 0 with very big file. 
This peak will be solved in next release and less memory will be needed.

6) Input Data:
-----------

A SAM file produced by an aligner with paired reads (BWA, Bowtie) and compliant with the SAM format. The reads must be paired.
If the input file is striped the buffer cache can be used.

7) Outputs: 
--------

Output are bam files per chromosome.
A bam files for discordant reads.
A bam file for unmapped reads.

Up to now the only way to uncompress the results: samtools view -Sh chrN.bam > chrN.sam <br />

8) MPI version:
------------

A MPI version should be installed first. We have tested the program with different MPI flavour: OpenMPI 1.10, MVAPICH.

9) Compiler: 
---------

A C compiler must be present. We have tested the programm with GCC and Intel Compiler. 

10) Sample test:
------------

We furnish a sample sam to sort for testing your installation.
We test it with from 1 to 8 jobs and 2 nodes with a normal network and file system. 

11) Lustre configuration:
-------------------

For Lustre and for our experiment the max_cached_mb is 48 Gb. That is 1/3 of the total RAM on each nodes for a node with 16 cores and 128 GB of memory
it makes 2.5 GB per jobs. 

According to Lustre documentation you can reach 3/4 of the total node memory. 

To see the amount you really cache use "top" on a node and look at the cached figures.

If you exceed the max_cached_mb the swap will enter in game and decrease performances. 
For the storage of the data we chose the largest striping of 2.5 GB.

At TGCC the striping factor is 128 (number of OSSs servers) and the striping size is maximum 2.5 GB. 

Before running and compile you must tell the programm how the input data is striped on Lustre and how to stripe the output data.
The parallel writing and reading are done via the MPI finfo structure. 

!!!We recommand to test different parameters before setting them once for all.

For reading and set the Lustre buffer size.
To do that edit the code of mpiSort.c and change the parameters in the header part. 
After tuning parameters recompile the application.

If you are familiar with MPI IO operation you can also test different commands collective, double buffered, data sieving.
In file write2.c in the function read_data_for_writing and writeSam, writeSam_unmapped, writeSam_discordant

The default parameters are for 128 OSS servers, with 2.5GB striping unit (maximum).
We do data sieving reading and a colective write is done with 128 servers.
The default are unharmed for other filesystem.

Tune this parameters according to your configuration.

12) File system management:
-----------------------

The cache optimization is critical for IO operations. 
The sorting algorithm does repeated reads to avoid disk access as explain in the Lustre documentation chapter 31.4. 
The cache keep the data of individual jobs in RAM and drastically improve IO performance.

You have the choice to use OSS cache or client (computer node) cache. Setting the cache at servers level is explained chapter 31.4.3.1.
Setting client cache is done via RPC as explain in chapter 31.4.1. Using the cache at client (computing node) side is recommended.

Here are parameters we set at TGCC for our experiments for client cache

Client side
 
max_cached_mb: 48301 <br />
max_read_ahead_mb=40 <br />
max_pages_per_rpc=256 <br />
max_rpcs_in_flight=32 <br />

The max_cached_mb tells how much cache can use the client. 

From our experiment on computing nodes with 16 cpu and 128 GB of RAM, a cache of 48 GB is enough <br />
This cache is 1/3 of the total memory on server and approximately 2.5 Gb per cpu. <br />
 
For NFS users there is a bug.
https://www.open-mpi.org/community/lists/users/2016/06/29434.php

a patch is available here:
https://trac.mpich.org/projects/mpich/attachment/ticket/2338/ADIOI_NFS_ReadStrided.patch

13) Results:
---------

To see results and speed-ups please check out those links

2015: This link was the first presentation of the tools 
http://devlog.cnrs.fr/_media/jdev2015/poster_jdev2015_institut_curie_hpc_sequencage_jarlier.pdf?id=jdev2015%3Aposters&cache=cache

2016: OpenSFS conference Lustre LAD 2016 
http://www.eofs.eu/_media/events/lad16/03_speedup_whole_genome_analysis_jarlier.pdf


14) Job rank:
----------

The number of jobs or jobs rank rely on the size of input data and the cache buffer size.
For our development on TGCC the cache size is 2.5GB per job.
 
If you divide the size input by the buffer cache you got the number of jobs you need. 
Of course the more cache you have the less jobs you need.

15) Example of code:
-----------------
How to launch the program with 

!/bin/bash
MSUB -r mpiSORT_HCC1187_20X_380cpu <br />
MSUB -@ frederic.jarlier@curie.fr:begin,end <br />
MSUB -n 380 <br />
MSUB -T 6000 <br />
MSUB -q large <br />
MSUB -o $SCRATCHDIR/ngs_data/ERROR/output_%I.o <br />
MSUB -e $SCRATCHDIR/ngs_data/ERROR/erreur_%I.e <br />
MSUB -w <br />

mpiSORT_BIN_DIR=$SCRATCHDIR/script_files/mpi_SORT <br />
BIN_NAME=psort <br />
FILE_TO_SORT=$SCRATCHDIR/ngs_data/HCC1187C_20X.sam <br />
OUTPUT_DIR=$SCRATCHDIR/ngs_data/sort_result/20X/ <br />
FILE=$SCRATCHDIR/ngs_data/sort_result/psort_time_380cpu_HCC1187_20X.txt <br />
mprun $mpiSORT_BIN_DIR/$BIN_NAME $FILE_TO_SORT $OUTPUT_DIR -q 0 <br />

16) Options 
----------

the -q option is for quality filtering. <br />
the -n for sorting by name <br />

17) Improvements
---------------

1) Optimize memory pressure on job rank 0 <br />
2) Manage single reads <br />
3) Mark or remove duplicates <br />
4) Make a pile up of the reads <br /> 
5) Write SAM files <br />
6) Propose an option to write a big SAM file <br />


