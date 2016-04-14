Requirements:
-------------

This programm is intended to run on supercomputer architectures equiped with Infiniband network and parallel 
filesystem such as Lustre, GPFS,...

The presented version of the program has been tested on France Genomic cluster of TGCC (Très Grand Centre de Calcul) of CEA (Bruyeres le Chatel, France). 

!!!!! DON'T RUN THIS PROGRAM ON NETWORK FILESYSTEM AKA NFS IT DOESN'T WORK!!!!

!!!!! DON'T RUN THIS PROGRAM IF YOU DON'T HAVE INFINIBAND OR LOW LATENCY NETWORK !!!!


You need to have a good understanding of your infrastructure before running and tuning this code.

Contact us if you need information.

Input Data:
----------

a SAM file produced by an aligner (BWA, Bowtie) and compliant with the SAM format. The reads should be paired.  

MPI version:
------------

A MPI version should be installed first. We have tested the program with different version: OpenMPI 1.10, MVAPICH.

Compiler: 
---------

A C compiler must be present also. We have tested the programm with GCC and Inter Compiler. 


Configuration:
--------------

The programm make an intensive use of reading buffer at the parallel filesystem level. 
According to the size of the Lustre or GPFS system the buffer size could vary. 
At TGCC the striping factor is 128 (number of OSSs servers) and the striping size is 2.5Gb. 
Before running and compile you must tell the programm how the input data is striped on Lustre and how to stripe the output data.

The reading and writing are done in different ways so the configuration varies between the two.

The parallel writing and reading are done via the MPI finfo structure. 

!!!We recommand to test different parameters before setting them once for all. 
Those parameter are independant of the file size you want to sort.    

1) for reading and set the Lustre buffer size.

To do that edit the code of parallelMergeSort.c and change the parameters from the line 169 to 182. 

 
2) for the writing part

Parameter for the writing part are located in the write2.c file from the line 2333 to 2340. 
Writing are done with collective operation so you have to tell how many buffer nodes you have.
After writing tell the programm to come back to parameters reading in the write2.c file from the line 2367 to 2373.

3) If you are familiar with MPI IO operation you can also test different commands collective, double buffered, data sieving.

In write2.c: line 2362 and 629.
 

Cache tricks and sizes:
----------------------

The cache optimization is critical for IO operations. 
The sorting algorithme does repeated reads to avoid disk access as explain in the Lustre documentation chapter 31.4. 
The cache keep the data of individual jobs in RAM and drastically improve IO performance.

You have the choice to use OSS cache or client cache. Setting the cache at servers level is explained chapter 31.4.3.1.   
Setting client cache is done via RPC as explain in chapter 31.4.1.

Here are parameters we set for our experiments for client cache

max_cached_mb: 48301
max_read_ahead_mb=40
max_pages_per_rpc=256
max_rpcs_in_flight=32

Form our experiments on Lustre with 64 to 128 OSSs the cache is 2GB with 16 OSSs the cache is 1GB.
Some parameters need root access ask your IT how to do this.

Job rank:
---------

The number of jobs or jobs rank rely on the size of input data and the cache buffer size.
For our development on TGCC the cache size is 2GB per job. 
If you divide the size input by the buffer cache you got the number of jobs you need. 
Of course the more cache you have the less jobs you need.

Default parameters:
-------------------

The default parameters are for 128 OSS servers, with 2GB striping unit (maximum).
We do data sieving reading and a colective write is done with 128 servers.

Tune this parameters according to your configuration.


Compilation:
------------

To compile the program modify the makefile to tell where the mpicc is located. Then type make in the installation folder.


Example of code:
-----------------
How to launch the program with 


!/bin/bash                                                                                                                                                                     
MSUB -r mpiSORT_HCC1187_20X_380cpu                                                                                                                                      
MSUB -@ frederic.jarlier@curie.fr:begin,end                                                                                                                                    
MSUB -n 380                                                                                                                                                                    
MSUB -T 6000                                                                                                                                                                   
MSUB -q large                                                                                                                                                                  
MSUB -o $SCRATCHDIR/ngs_data/ERROR/output_%I.o                                                                                                        
MSUB -e $SCRATCHDIR/ngs_data/ERROR/erreur_%I.e                                                                                                        
MSUB -w                                                                                                                                                                        
set -x
cd ${BRIDGE_MSUB_PWD}

mpiSORT_BIN_DIR=$SCRATCHDIR/script_files/mpi_SORT
BIN_NAME=psort
FILE_TO_SORT=$SCRATCHDIR/ngs_data/HCC1187C_20X.sam
OUTPUT_DIR=$SCRATCHDIR/ngs_data/sort_result/20X/
FILE=$SCRATCHDIR/ngs_data/sort_result/psort_time_380cpu_HCC1187_20X.txt
ccc_mprun $mpiSORT_BIN_DIR/$BIN_NAME $FILE_TO_SORT $OUTPUT_DIR -q 0

Options 
-------

the -q option is for quality filtering.


Future developments
-------------------

1) Add options to manage cache from command lines
2) Mark and remove duplicates
3) Make a pile up of the reads 



