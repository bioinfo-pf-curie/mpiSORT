/*
   mpiSORT
   Copyright (C) 2016-2017 Institut Curie, 26 rue d'Ulm, Paris, France

   mpiSORT is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   mpiSORT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser Public License
   along with mpiSORT.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
   Module:
     bufferized_read.c

   Authors:
 	Frederic Jarlier from Institut Curie
	Nicolas Joly from Institut Pasteur
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/

#ifndef WRITE_H
#define WRITE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <mpi.h>

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

void writeSam(
		int rank,
		char* output_dir,
		char* header,
		size_t local_readNum,
		size_t total_num_read,
		char* chrName,
		Read* chr,
		int total_num_proc,  //this is the number of proc in split communication
		MPI_Comm split_comm,
		int master_rank,
		char *file_name,
		MPI_File in,
		MPI_Info finfo,
		int compression_level,
		size_t* new_offset_dest,
		size_t* new_offset_source,
		int* new_read_size,
		int* new_rank,
		char* data,
		size_t offset_data_in_file);


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
		int num_proc, MPI_Comm split_comm, char *file_name, MPI_File in, MPI_Info finfo, int compression_level, char *data,
		size_t offset_data_in_file, int write_sam);

#endif
