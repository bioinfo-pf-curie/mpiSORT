# Advanced tuning of the source code for parallel filesystems

As the programs use MPI functions for reading and writing you can take advantage of a parallel file system. To speed-up reading and writing you can set the striping of your data. The striping tells the number of file servers you want to use and the size of data blocks on each server. You can compare the striping of the file with the mapping process used in Hadoop. This is the way your data are distributed among the servers of your paralell filesystem. This kind of optimization accelerate drastically the IO operations.

Standard software for sorting SAM file don't take into account the underlying distributed filesystem and low-latency interconnexion while `mpiSORT` does. Indeed, MPI improves IOs by means of parallelization. On parallel filesystem such as [Lustre](http://lustre.org/) or [GPFS](https://en.wikipedia.org/wiki/IBM_Spectrum_Scale), one way of accelerating IOs is to stripe them across servers. You can choose the number of servers with the striping factor and the size of chunks with the striping size. 

## Lustre configuration

Before compiling and running the program, you must specify how the input data are striped on [Lustre](http://lustre.org/) and how to stripe the output data. The parallel writing and reading are done via the MPI finfo structure. In order to modify the code according to your filesystem configuration performs as follows.


For the reading part there are two blocks of code in the `mpiSort.c` in which you can set the [Lustre](http://lustre.org/) buffer size, these blocks of code being surrounded by the following tags:

```
// BEGIN> FINE TUNING FINFO FOR WRITING OPERATIONS

// END> FINE TUNING FINFO FOR WRITING OPERATIONS


For writing part there are two blocks of codes in the `write.c` in which you can set the [Lustre](http://lustre.org/) buffer size, these blocks of code being surrounded by the following tags:

```
// BEGIN> FINE TUNING FINFO FOR WRITING OPERATIONS

// END> FINE TUNING FINFO FOR WRITING OPERATIONS
```

These parameters are harmless for other filesystem.

Note that, if you are familiar with MPI IO operations you can also test the different commands such as *collective*, *double buffered*, *data sieving*. they can be changed in the following functions `writeSam`, `writeSam_unmapped` and `writeSam_discordant` from `write.c`.

 At [TGCC France Génomique](https://www.france-genomique.org/plateformes-et-equipements/plateforme-tgcc-arpajon/) (Broadwell architecture), with a striping factor of 128 (number of OSSs servers) and the a striping size of 2.5 GB, a SAM file of 1.2TB could be loaded in memory in less than 1 minute.

We strongly recommand to test different parameters to assess what are the optimal setting for your configuration.
