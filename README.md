Objective
---------

Sorting big NGS data file Version 1.0. 

Sections:
-------

1) Release notes <br />
2) Installation <br />
3) Algorithm <br />
4) Architectures <br />
5) Memory usage <br />
6) Inputs <br />
7) Outputs <br />
8) Compiler <br />
9) Sample test  <br />
10) Configuration <br />
11) Parallel file system configuration <br />
12) Network file system configuration <br />
13) Results <br />
14) Example of code <br />
15) Options <br />
16) Improvements and future work <br />

==============================================================================


1) Release notes 
-------------

Release 1.0 from 11/10/2016

1) In this release we use the local data buffer for re-reading instead of the linux client cache. <br />
Because for slow network the time to get the data back from the client or server cache can be long. <br />
Now re-reading the data is intantaneous. <br />
2) Add documentation and free some variables. <br />

Release 1.0 from 06/10/2016

1) The previous version didn't sort the offset destination before the shuffle. This bug is fixed. <br />
2) New packaging with autools for distribution. <br />

Release 0.9 from 29/07/2016

1) Bug fix when data sample is too small and some jobs have nothing to sort <br />
2) Bug fix when jobs under 6 <br />
3) We add new option -n for sorting by name <br />
4) the parallel merge sort has been replaced with bitonic sorting (25% percent gain in speed-up) <br />
5) Refactoring and cleaning of the code <br />

2) Installation
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

As the programs use MPI fonctions for reading and writing you can take advantage of a parallel file system. To speed-up reading and writing you can set the striping of your data the striping tells the number of file servers you want to use and the size of data blocks on each server. You can compare the striping of the file with the mapping process use in Hadoop. This is the way your data are distributed among the servers of your file system. This kind of optimizations accelerate drastically the IO operations.<br />

Ordinary softwares (Samtools, Sambamba, Picard,... ) doesn't take into account the underlying distributed file system and low latency interconnexion when MPI does. <br />

4) Architectures:
-------------

We have tested the programm on different architectures.

In the Institut Curie we have a NFS and a 10Gb network. The France Genomic cluster of TGCC (Très Grand Centre de Calcul) of CEA (Bruyeres le Chatel, France) is equiped with Lustre. <br />

Because of the development design the programm is optimized for HPC architecture. This programm runs better on low latency network and parallel file system. <br />

Contact us if you need information.

5) Memory
---------

The total memory used during the sorting is around three times the size of the SAM file. <br />

6) Input Data:
-----------

A SAM file produced by an aligner with paired reads (BWA, Bowtie) and compliant with the SAM format. The reads must be paired. <br />

7) Outputs: 
--------

Output are bam files per chromosome, a bam file for discordant reads and a bam file for unmapped reads. <br />

To index the bam file use tabix like this: <br />
tabix -p sam chrN.bam <br />
(samtools index doesn't work yet)<br />

To uncompress the results:<br />
samtools view -Sh chrN.bam > chrN.sam <br />

Or rename the file with suffix .gz and 
bgzip -d chrN.gz > chrN.sam

8) MPI version:
------------

A MPI version should be installed first. We have tested the program with different MPI flavour: OpenMPI 1.10.0 and 1.8.3.

9) Compiler: 
---------

A C compiler must be present. We have tested the programm with GCC and Intel Compiler. 

10) Sample test:
------------

We furnish a sample sam to sort for testing your installation.
We test it with from 1 to 8 jobs and 2 nodes with a normal network and file system. 

11) Parallel file system configuration:
-----------------------------

MPI improves IOs by means of parallelization. 

On parallel file system like Lustre, GPFS one way of accelerating IO is to stripe them across servers. 
You chose the number of servers with the striping factor and the size of chunks with the striping size. 

Before running and compile you must tell the programm how the input data is striped on Lustre and how to stripe the output data.
The parallel writing and reading are done via the MPI finfo structure in the first line of mpiSort.c.

!!!We recommand to test different parameters before setting them once for all.

For reading and set the Lustre buffer size.
To do that edit the code of mpiSort.c and change the parameters in the header part. 
After tuning parameters recompile the application.

If you are familiar with MPI IO operation you can also test different commands collective, double buffered, data sieving.
In file write.c in the function writeSam, writeSam_unmapped, writeSam_discordant

The default parameters of the programm are unharmed for other filesystem.

At TGCC the striping factor is 128 (number of OSSs servers) and the striping size is maximum 2.5 GB a SAM file of 1.2TB could be loaded in memory in less than 1 minute.

12) Network File system configuration:
-------------------------------

We don't recommend to use MPI with NFS (it works but it does not scale very well). 

For NFS users there is a bug in openMPI.
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

14) Example of code:
-----------------
How to launch the program with Torque:

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

15) Options 
----------

the -q option is for quality filtering. <br />
the -n for sorting by name <br />

16) Improvements
---------------

1) Generate index with the bam 
2) Optimize memory pressure on job rank 0 (do the indexing in the bitonic sort)<br />
3) Manage single reads <br />
4) Mark or remove duplicates <br />
5) Make a pile up of the reads <br /> 
6) Write SAM files per chromosom<br />
7) Propose an option to write a big SAM/BAM file <br />


