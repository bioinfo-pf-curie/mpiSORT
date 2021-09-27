/*
   This file is part of mpiSORT
   
   Copyright Institut Curie 2020
   
   This software is a computer program whose purpose is to sort SAM file.
   
   You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
   
   The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
   
   The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/

/*
   Module:
     writeUtils.h

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/

#ifndef WRITE_UTILS_H
#define WRITE_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parser.h"

extern size_t g_wu_readNum;
extern size_t g_wu_totalReadNum;
extern size_t g_wu_master;
extern MPI_Comm COMM_WORLD;
extern MPI_Status status;

size_t file_get_size(char * path);


static inline void unpack_buffer(void * dest, void * src, int count, void **indices, int *blocklens, size_t elem_size)
{
	int i;

	size_t curr_off = 0;
	for( i = 0; i < count; i++)
	{
		memcpy(dest + (size_t)indices[i], src + curr_off, blocklens[i] * elem_size);
		curr_off += blocklens[i] * elem_size;
	}

}


void create_read_dt(int rank, int num_proc, int *ranks, int* buffs, char** data, MPI_Datatype* dt, size_t readNum);

void * create_send_datatype_for_size(int rank, int size, size_t *num_reads_by_procs, int **dest_size,
		int k, size_t* pack_size, int** recv_index);

void * create_send_datatype_for_offsets(int rank, int size, size_t *num_reads_by_procs, size_t **dest_offsets,
		int k, size_t* packed_size, int** recv_index);

void * create_send_datatype_for_reads(int rank, int size, size_t *buffs, char** data, int k, size_t* packed_size, int** recv_index);

size_t get_send_size(int rank, int size, size_t* buffs, size_t** send_size, int count, int k);




void send_size_t_all_to_master(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs,
		size_t *start_size_per_job, size_t *all_data, size_t *data);
void send_size_t_all_to_master_bitonic(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs,
		size_t *start_size_per_job, size_t *all_data, size_t *data, size_t start_index);
void send_size_t_all_to_master_bitonic_V2(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs,
		size_t *start_size_per_job, size_t *all_data, size_t *data, size_t start_index);
void send_size_t_master_to_all(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs, size_t *start_size_per_job, size_t *all_data, size_t *data);
void send_int_all_to_master(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs, size_t *start_size_per_job, int *all_data, int *data);
void send_int_master_to_all(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs, size_t *start_size_per_job, int *all_data, int *data);


int badCount(int k, int n);

//Maybe useful in some other file
void assert_sorted_size_t(size_t * buffer, size_t number);

#endif
