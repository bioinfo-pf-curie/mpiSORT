#ifndef WRITE_UTILS_H
#define WRITE_UTILS_H

#include <stdio.h>
#include <stdlib.h>

#include "parser.h"

extern size_t g_wu_readNum;
extern size_t g_wu_totalReadNum;
extern size_t g_wu_master;
extern MPI_Comm COMM_WORLD;
extern MPI_Status status;

void create_read_dt(int rank, int num_proc, int *ranks, int* buffs, char** data, MPI_Datatype* dt, size_t readNum);

int create_send_datatype_for_size(int rank, int size, size_t *num_reads_by_procs, int **dest_size,
		int k, MPI_Datatype* dt, int** recv_index);

int create_send_datatype_for_offsets(int rank, int size, size_t *num_reads_by_procs, size_t **dest_offsets,
		int k, MPI_Datatype* dt, int** recv_index);

int create_send_datatype_for_reads(int rank, int size, size_t *buffs, char** data, int k, MPI_Datatype* dt, int** recv_index);

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
