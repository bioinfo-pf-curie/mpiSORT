
Release 1.0 from the 09/01/2020  <br />

1) remove character "\0" during parse paired

Release 1.0 from the 31/03/2019  <br />

1) fix a corner case when reads distribution is unbalanced.  <br />

Release 1.0 from the 28/06/2018

1) Creation of a google group
	https://groups.google.com/forum/#!forum/hpc-bioinformatics-pipeline <br />

2) For improvement on Lustre file system removed the “nosuid” mount option. <br />


Release 1.0 from 30/04/2018

1) fix a memory leak with khash in write.c <br />

Release 1.0 from 16/02/2018

1) fix a bug in the offset of the first read. <br />
Change the argument of init_goff function replace headerSize with hsiz in mpiSort.c.<br />

2) next steps are: <br />
    -Investigate why Valgrind complains about MPI_mem_alloc in bitonic sort part?<br />
    -Marking duplicates. <br />

Release 1.0  from 30/08/2017

1) fix the condition for quality treshold. (parser.c line 179)
Now we include read with mapping quality equal to zero (when no quality is specified). <br />

Release 1.0  from 14/04/2017

1) Add citations.  <br />
2) Improve the memory pressure (see memory usage).  <br />
3) Remove some assert and add comments.  <br />
4) Malloc_trim is no op on BSD or OSX.  <br />
5) First tests (see results).  <br />

Release 1.0 from 28/03/2017

1) We have redesign the communications during the sorting. <br />
Communication are optimized for power of 2 in processor number.  <br />
This has been done to scale gigabit ethernet network and for a better load balancing of the computing and memory pressure. <br />
With processors in power of 2 we avoid extra communication. <br />
All Results are the same than previous release. <br />

New tests and benchmarking will follow.

2) Remove memory leaks. <br />
Comments added. <br />

Release 1.0 from 25/11/2016

1) Bug fix at the beginning of the programm during the parsing of the SAM.<br />
2) Need to test if the MPI sctructures scale very large buffer (above 5Gb per job).<br />

Release 1.0 from 23/11/2016

1) Remove the limit in  MPI_file_read_at  of 1 Gb per read buffer.<br />
SAM data are first stored in intermediate buffer and copy to the main SAM buffer see in mpiSORT.c (line 344)<br />
2) Replace MPI_Create_struct with a local copy (in write.c line 2712) <br />
3) Cleaning of the code and add comments <br />
4) Tested with a beegfs scratch (4 MDS) and 10gb ethernet, it scales to 200Gb sam file and 80 cpu and takes 10mn for reading and 10mn sorting). <br />

Release 1.0 from 18/11/2016

1) Trimming of the memory after big or multiple free() with malloc_trim().<br /> 
Efficient with Linux but not tested on BSD of IOX. <br />
2) Remove memory leaks.<br />
3) Add/remove  comments.<br />

Release 1.0 from 04/11/2016

1) Due to MPI packing overhead we replace with a local copy (write.c). See in future if MPI_Unpack is stable.  <br />
2) File extension becomes  .gz (not bam) <br />
3) Discordant and unmapped are computed first.  <br />
4) Cleaning of memory.  <br />

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
