#ifndef WRITE_H
#define WRITE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "mpi.h"
#include "parser.h"
#include "assert.h"
#include "preWrite.h"
#include "diffuse.h"
#include "math.h"
#include "sys/mman.h"
#include "mpi_globals.h"
#include "qksort.h"
#include "write_utils.h"
#include "parabitonicsort.h"

/*
void writeSam(
		char* output_dir, char *file_name, MPI_File in, MPI_Info finfo,
		char* header,
		char* chrName, Read* reads, size_t* offsets, size_t readNum,
		int compression_level,
		MPI_Comm comm);
*/

void writeSam(size_t total_num_reads,
		int master_rank,
		int split_comm_size,
		int dimensions,
		char* output_dir,
		char *file_name,
		MPI_File in,
		MPI_Info finfo,
		char* header,
		char* chrName, Read* reads,
		size_t *local_dest_offsets_sorted,
		size_t *local_source_offsets_sorted,
		int *local_read_size_sorted,
		int *local_rank_sorted,
		size_t readNum,
		int compression_level,
		MPI_Comm comm);

size_t init_offset_and_size_free_chr(size_t* offset, int* size, Read* data_chr, int local_readNum);
void read_data_for_writing(int rank, int num_proc, size_t local_readNum, char *file_name,
		size_t *number_of_reads_by_procs, size_t *buffs_by_procs, char *** data,
		int *new_rank, int *new_size, size_t *new_offset, MPI_File in, MPI_Info finfo, MPI_Comm COMM_WORLD);

void bruckWrite(int rank, int num_proc,
		size_t local_readNum, size_t* number_of_reads_by_procs, int *new_rank,
		size_t *buffs_by_procs, char*** data2,
		size_t *new_offset, size_t*** data_offsets,
		int *new_size, int ***data_size
	);
void bruck_reads(int rank, int num_proc, size_t * buffs_by_procs, char** data2);
void bruck_offsets(int rank, int num_proc, int local_readNum, size_t* number_of_reads_by_procs, size_t ** data_offsets, int *new_rank, size_t* new_offset);
void bruck_size(int rank, int num_proc, size_t local_readNum, size_t* number_of_reads_by_procs, int ** data_size, int *new_rank, int *new_size);

void writeSam_discordant_and_unmapped(int split_rank, char* output_dir, char* header, size_t local_readNum, char* chrName, Read* chr,
		int num_proc, MPI_Comm split_comm, char *file_name, MPI_File in, MPI_Info finfo, int compression_level);

#endif
