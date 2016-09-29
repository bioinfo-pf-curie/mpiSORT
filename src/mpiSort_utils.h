
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

void get_coordinates_and_offset_source_and_size_and_free_reads(int rank, int *local_read_rank, size_t *coordinates,
		size_t* offset, int* size, Read* data_chr, int local_readNum);

size_t init_coordinates_and_size(int rank, int *local_reads_rank, size_t *local_reads_index,
		size_t* coordinates, int* size, Read* data_chr, int local_readNum);


void chosen_split_rank_gather_size_t(MPI_Comm split_comm, int rank, int num_proc, int master, size_t size, size_t *size_per_jobs,
		size_t *start_size_per_job, size_t *all_data, size_t *data, size_t start_index);
