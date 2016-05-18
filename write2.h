/*
   mpiSORT
   Copyright (C) 2016 Institut Curie, 26 rue d'Ulm, Paris, France

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
     write2.h

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/



typedef struct {
		//int  read_lenght;
		char *read;
	} Read_to_write;

size_t init_offset_and_size_free_chr(size_t* offset, int* size, Read* data_chr, int local_readNum);
void bruck_reads(int rank, int num_proc, size_t * buffs_by_procs, char** data2);
void bruck_offsets(int rank, int num_proc, int local_readNum, size_t* number_of_reads_by_procs, size_t ** data_offsets, int *new_rank, size_t* new_offset);
void bruck_size(int rank, int num_proc, size_t local_readNum, size_t* number_of_reads_by_procs, int ** data_size, int *new_rank, int *new_size);

void bruckWrite(int rank, int num_proc,
		size_t local_readNum, size_t* number_of_reads_by_procs, int *new_rank,
		size_t *buffs_by_procs, char*** data2,
		size_t *new_offset, size_t*** data_offsets,
		int *new_size, int ***data_size
	);

void read_data_for_writing(int rank, int num_proc, size_t local_readNum,
		char *file_name, size_t *number_of_reads_by_procs, size_t *buffs_by_procs,
		char *** data, int *new_rank, int *new_size, size_t *new_offset, MPI_File in, MPI_Info finfo, MPI_Comm COMM_WORLD);

void send_size_t_all_to_master(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs, size_t *start_size_per_job, size_t *all_data, size_t *data);
void send_int_all_to_master(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs, size_t *start_size_per_job, int *all_data, int *data);
void send_size_t_master_to_all(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs, size_t *start_size_per_job, size_t *all_data, size_t *data);
void send_int_master_to_all(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs, size_t *start_size_per_job, int *all_data, int *data);

void writeSam(int rank, char* output_dir, char* header, size_t readNum,
			char* chrName, Read* chr, size_t* offsets, int num_proc, MPI_Comm NEW_COMM_WORLD,
			char *file_name, MPI_File in, MPI_Info finfo, int compression_level);

void writeSam_unmapped(int rank, char* output_dir, char* header, size_t readNum,
			char* chrName, Read* chr, int num_proc,
				MPI_Comm NEW_COMM_WORLD, char *file_name,
					MPI_File in, MPI_Info finfo, int compression_level);

void writeSam_discordant(int rank, char* output_dir, char* header, size_t readNum,
			char* chrName, Read* chr, int num_proc,
				MPI_Comm NEW_COMM_WORLD, char *file_name,
					MPI_File in, MPI_Info finfo, int compression_level);

int issort(void *data, size_t size, size_t esize, int (*compare)(const void *key, const void *key2));
int qksort(void *data, size_t size, size_t esize, size_t i, size_t k, int (*compare)(const void *key1, const void *key2));
int badCount(int k, int n);


/*
 * Definitions for the modified Bruck
 * Part in write2.c
 */
void create_DT(int rank, int numprocs, int k, MPI_Datatype oldtype, MPI_Datatype* newtype);
int get_block_to_read_size(int rank, int size, char* file_name, int* buffsize, int* ranks, int * buffs);
void create_read_dt(int rank, int size, int *ranks, int *buffs, char** data, MPI_Datatype* dt, size_t readNum);

int create_send_datatype_for_reads(int rank, int size, size_t *buffs, char** data, int k, MPI_Datatype* dt, int** recv_index);

size_t get_send_size(int rank, int size, size_t* buffs, size_t** send_size, int count, int k);
int badCount(int k, int n);

int create_send_datatype_for_offsets(int rank, int size, size_t *num_reads_by_procs, size_t **dest_offsets,
		int k, MPI_Datatype* dt, int** recv_index);


int create_send_datatype_for_size(int rank, int size, size_t *num_reads_by_procs, int **dest_offsets,
		int k, MPI_Datatype* dt, int** recv_index);

void ParallelBitonicSort(int my_rank, int dimension, size_t *local_list, size_t *local_index, size_t list_size, size_t zero_padding);

void send_size_t_all_to_master_bitonic(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs,
		size_t *start_size_per_job, size_t *all_data, size_t *data, size_t start_index);

void send_size_t_all_to_master_bitonic_V2(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs,
		size_t *start_size_per_job, size_t *all_data, size_t *data, size_t start_index);

