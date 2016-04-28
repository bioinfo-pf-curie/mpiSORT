Requirements:
-------------

This programm is intended to run on supercomputer architectures equiped with Infiniband network and parallel 
filesystem such as Lustre, GPFS,...

The presented version of the program has been tested on France Genomic cluster of TGCC (Très Grand Centre de Calcul) of CEA (Bruyeres le Chatel, France). 

!!!!! THIS PROGRAM NEEDS A LOW LATENCY NETWORK !!!! <br />
!!!!! THIS PROGRAM NEEDS A PARALLEL FILE SYSTEM !!!!

Optimizations are everywhere :)

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

The programm make an intensive of cache buffer at the client level. The cache size depends of the variable you give to Lustre.
In our experiement the max_cached_mb is 48 Gb. That is 1/3 of the total RAM on each nodes for a node with 16 cores and 128 GB of memory
it makes 2.5 GB per jobs. 

According to Lustre documentation you can reach 3/4 of the total node memory. 

To see the amount you really cache use "top" on a node and look at the cached figures.

If you exceed the max_cached_mb the swap will enter in game and decrease performances. 
 
For the storage of the data we chose the largest striping of 2 GB.

At TGCC the striping factor is 128 (number of OSSs servers) and the striping size is 2 GB. 

Before running and compile you must tell the programm how the input data is striped on Lustre and how to stripe the output data.

The reading and writing are done in different ways so the configuration varies between the two.

The parallel writing and reading are done via the MPI finfo structure. 

!!!We recommand to test different parameters before setting them once for all.    

1) for reading and set the Lustre buffer size.

To do that edit the code of parallelMergeSort.c and change the parameters from the line 169 to 182. 

MPI_Info_create(&finfo);<br />
MPI_Info_set(finfo,"striping_factor","128"); OSS number <br />
MPI_Info_set(finfo,"striping_unit","2684354560"); 2G striping <br />
MPI_Info_set(finfo,"ind_rd_buffer_size","2684354560"); 2gb buffer <br />
MPI_Info_set(finfo,"romio_ds_read","enable"); for data sieving reading <br />
		
MPI_Info_set(finfo,"nb_proc","128"); for collective buffering write <br />
MPI_Info_set(finfo,"cb_nodes","128"); <br />
MPI_Info_set(finfo,"cb_block_size","2147483648"); <br /> 
MPI_Info_set(finfo,"cb_buffer_size","2147483648"); <br />

After tuning parameters recompile the application.
 
2) for the writing part

Parameter for the writing part are located in the write2.c file from the line 2333 to 2340. 

MPI_Info_set(finfo,"striping_factor","128"); // 128 OSS Lustre <br />
MPI_Info_set(finfo,"striping_unit","268435456"); //256 MBytes striping <br />

MPI_Info_set(finfo,"nb_proc","64"); // MPI_write <br />
MPI_Info_set(finfo,"cb_nodes","64"); // MPI_write buffer nodes<br />
MPI_Info_set(finfo,"cb_block_size","268435456"); /* 256 MBytes - should match FS block size */ <br />
MPI_Info_set(finfo,"cb_buffer_size","536870912"); /* 512 MBytes or multiple of cb_block_size(Optional) */ <br />

Writing are done with collective operation so you have to tell how many buffer nodes you have.


Then after writing tell the programm to come back to parameters reading in the write2.c file from the line 2367 to 2373.

MPI_Info_set(finfo,"striping_factor","128"); <br />
MPI_Info_set(finfo,"striping_unit","2684354560"); <br />

MPI_Info_set(finfo,"nb_proc","128"); <br />
MPI_Info_set(finfo,"cb_nodes","128"); <br />
MPI_Info_set(finfo,"cb_block_size","2684354560"); <br /> 
MPI_Info_set(finfo,"cb_buffer_size","2684354560"); <br />


Rmq: In future release those options will pe passed in argument. 

3) If you are familiar with MPI IO operation you can also test different commands collective, double buffered, data sieving.

In write2.c: line 2362 and 629.
 

Cache tricks and sizes:
----------------------

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
 

Job rank:
---------

The number of jobs or jobs rank rely on the size of input data and the cache buffer size.
For our development on TGCC the cache size is 2.5GB per job.
 
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

mpiSORT_BIN_DIR=$SCRATCHDIR/script_files/mpi_SORT <br />
BIN_NAME=psort <br />
FILE_TO_SORT=$SCRATCHDIR/ngs_data/HCC1187C_20X.sam <br />
OUTPUT_DIR=$SCRATCHDIR/ngs_data/sort_result/20X/ <br />
FILE=$SCRATCHDIR/ngs_data/sort_result/psort_time_380cpu_HCC1187_20X.txt <br />
ccc_mprun $mpiSORT_BIN_DIR/$BIN_NAME $FILE_TO_SORT $OUTPUT_DIR -q 0 <br />

Options 
-------

the -q option is for quality filtering.


Future developments
-------------------

1) Add options to manage cache optimizations from command lines <br />
2) Mark and remove duplicates <br />
3) Make a pile up of the reads <br /> 



