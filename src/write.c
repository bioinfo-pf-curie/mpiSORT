/*
   mpiSORT
   Copyright (C) 2016-2017 Institut Curie / Institut Pasteur

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
    write.c

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <malloc.h>

#include "compat.h"

#include "write.h"
#include "bgzf.h"
#include "bgzf.c"
#include "mpiSort_utils.h"
#include "parabitonicsort2.h"

size_t init_offset_and_size_free_chr(size_t* offset, int* size, Read* data_chr, int local_readNum)
{
	size_t dataSize = 0;
	int j;
	Read* chr = data_chr;
	Read* to_free = NULL;

	//we initialize offset source and size_source
	for(j = 0; j < local_readNum; j++){
		size[j] = 0;
		offset[j] = 0;
	}

	//first we are going to read in the source file
	for(j = 0; j < local_readNum; j++){
		//offset is the read size
		offset[j] = chr->offset_source_file;
		size[j] = (int)chr->offset; //read length
		dataSize += chr->offset;

		to_free = chr;
		chr = chr->next;

		free(to_free);
	}

	return dataSize;
}

//BRUCK FUNC
void bruckWrite(int rank, int num_proc,
		size_t local_readNum, size_t* number_of_reads_by_procs, int *new_rank,
		size_t *buffs_by_procs, char*** data2,
		size_t *new_offset, size_t*** data_offsets,
		int *new_size, int ***data_size)
{

	bruck_reads(rank, num_proc, buffs_by_procs, *data2);
	bruck_size(rank, num_proc, local_readNum, number_of_reads_by_procs, *data_size, new_rank, new_size);
	bruck_offsets(rank, num_proc, local_readNum, number_of_reads_by_procs, *data_offsets, new_rank, new_offset);

}

void bruckWrite2(
		int rank,
		int num_proc,
		size_t local_readNum,
		size_t* number_of_reads_by_procs,
		int *new_rank,
		size_t *new_dest_offset,
		size_t*** data_dest_offsets,
		int *new_size,
		int ***data_size,
		size_t *new_source_offset,
		size_t*** data_source_offsets,
		int *new_dest_rank,
		int*** dest_rank
		){

	//we trade the size
	bruck_size(rank, num_proc, local_readNum, number_of_reads_by_procs, *data_size, new_rank, new_size);
	bruck_offsets(rank, num_proc, local_readNum, number_of_reads_by_procs, *data_dest_offsets, new_rank, new_dest_offset);
	bruck_offsets(rank, num_proc, local_readNum, number_of_reads_by_procs, *data_source_offsets, new_rank, new_source_offset);
	bruck_size(rank, num_proc, local_readNum, number_of_reads_by_procs, *dest_rank, new_rank, new_dest_rank);
}

void bruckWrite3(
		int rank,
		int num_proc,
		size_t local_readNum,
		size_t* number_of_reads_by_procs,
		int *new_rank,
		size_t *new_source_offset,
		size_t ***data_source_offsets,
		int *new_dest_rank,
		int ***dest_rank,
		size_t *new_reads_coordinates,
		size_t ***read_coordinates,
		int *new_reads_size,
		int ***read_size,
		int *new_source_rank,
		int ***source_rank,
		size_t *new_dest_offset,
		size_t ***dest_offset
		){


	bruck_offsets(rank, num_proc, local_readNum, number_of_reads_by_procs, *read_coordinates, new_rank, new_reads_coordinates);
	bruck_offsets(rank, num_proc, local_readNum, number_of_reads_by_procs, *dest_offset, new_rank, new_dest_offset);
	bruck_offsets(rank, num_proc, local_readNum, number_of_reads_by_procs, *data_source_offsets, new_rank, new_source_offset);
	bruck_size(rank, num_proc, local_readNum, number_of_reads_by_procs, *dest_rank, new_rank, new_dest_rank);
	bruck_size(rank, num_proc, local_readNum, number_of_reads_by_procs, *source_rank, new_rank, new_source_rank);
	bruck_size(rank, num_proc, local_readNum, number_of_reads_by_procs, *read_size, new_rank, new_reads_size);
}

void bruck_reads(int rank, int num_proc, size_t * buffs_by_procs, char** data2)
{
	MPI_Comm comm = COMM_WORLD;

	int k, m, srank, rrank;
	MPI_Datatype dt_send, dt_recv;
	size_t *recv_size_by_proc	= malloc(sizeof(size_t));
	size_t *send_size_by_proc	= malloc(sizeof(size_t));
	int *recv_index				= malloc(sizeof(int));
	char* interbuff 			= malloc(1);
	size_t total, send_total;
	int packsize;
	double time;

	int count;

	time = MPI_Wtime();
	for(k=1; k<num_proc; k<<=1)
	{
		srank = (rank - k + num_proc) % num_proc;	//Rank to send to
		rrank = (rank + k) % num_proc;	//Rank to recv from

		int count = badCount(k, num_proc);
		recv_index = realloc(recv_index, sizeof(int) * count);

		count = create_send_datatype_for_reads(rank, num_proc, buffs_by_procs, data2, k, &dt_send, &recv_index);
		MPI_Type_commit(&dt_send);

		send_size_by_proc = realloc(send_size_by_proc, count*sizeof(size_t));
		recv_size_by_proc = realloc(recv_size_by_proc, count*sizeof(size_t));

		send_total = get_send_size(rank, num_proc, buffs_by_procs, &send_size_by_proc, count, k);
		MPI_Pack_size(1, dt_send, comm, &packsize);

		assert(packsize == send_total);

		MPI_Sendrecv(send_size_by_proc, count, MPI_LONG_LONG_INT, srank, 0,
				recv_size_by_proc, count, MPI_LONG_LONG_INT,
				rrank, 0, comm, MPI_STATUS_IGNORE);

		total = 0;
		for(m = 0; m < count; m++)
		{
			total += recv_size_by_proc[m];
		}

		interbuff = realloc(interbuff, total+1);
		memset(interbuff,0,(total+1)*sizeof(char));

		MPI_Sendrecv(MPI_BOTTOM, 1, dt_send, srank, 0,
				interbuff, total, MPI_PACKED, rrank, 0, comm, MPI_STATUS_IGNORE);

		for ( m = 0; m < count; m++){
			// according to the recieve size
			if (data2[recv_index[m]]){
				data2[recv_index[m]] = realloc(data2[recv_index[m]], recv_size_by_proc[m]+1);
				data2[recv_index[m]][recv_size_by_proc[m]] = 0;
			}
		}

		MPI_Aint indices[count];
		int blocklens[count];
		MPI_Datatype oldtypes[count];

		for (m = 0; m < count; m++){

			blocklens[m] = (int)recv_size_by_proc[m];
			MPI_Get_address(data2[recv_index[m]], &indices[m]);
			oldtypes[m] = MPI_CHAR;
		}

		//Create structure of recieve type
		MPI_Type_create_struct(count, blocklens, indices, oldtypes, &dt_recv);
		MPI_Type_commit(&dt_recv);

		int pos=0;
		MPI_Unpack(interbuff, total, &pos, MPI_BOTTOM, 1, dt_recv, comm);

		for(m = 0; m<count; m++)
		{
			buffs_by_procs[recv_index[m]] = strlen(data2[recv_index[m]])+1;
		}
		MPI_Barrier(comm);
		MPI_Type_free(&dt_recv);
		MPI_Type_free(&dt_send);
		count = 0;
	}

	free(interbuff);
	free(recv_index);
	free(recv_size_by_proc);
	free(send_size_by_proc);
}

void bruck_offsets(int rank, int num_proc, int local_readNum, size_t* number_of_reads_by_procs, size_t ** data_offsets, int *new_rank, size_t* new_offset)
{
	MPI_Comm comm = COMM_WORLD;

	int k, m, j, srank, rrank;
	MPI_Datatype dt_send;

	size_t *recv_size_by_proc=NULL, *send_size_by_proc=NULL;
	int *recv_index=NULL;
	size_t total, send_total;
	int packsize;
	double time;

	int count;

	time = MPI_Wtime();

	for (m = 0; m < num_proc; m++){
		number_of_reads_by_procs[m] = 0;

	}

	for(m = 0; m < local_readNum; m++){
		number_of_reads_by_procs[new_rank[m]]++;
	}

	size_t **data_offsets2 = malloc(sizeof(size_t *)*num_proc);
	//we initialize data_offsets
	for(m = 0; m < num_proc; m++){
		data_offsets[m] = NULL;
		data_offsets2[m] = calloc( number_of_reads_by_procs[m], sizeof(size_t));
	}

	size_t *read_by_proc_offset = calloc(num_proc, sizeof(size_t));

	//we give values data_offsets
	for(m = 0; m < local_readNum; m++){

		// Phase one of bruck shift of the (rank+i)%size and (rank-i+size)%size
		data_offsets2[new_rank[m]][read_by_proc_offset[new_rank[m]]] = new_offset[m];
		read_by_proc_offset[new_rank[m]]++;
	}

	for (j = 0; j < num_proc; j++){
		data_offsets[(rank+j)%num_proc] = data_offsets2[(rank-j+num_proc)%num_proc];
		number_of_reads_by_procs[(rank+j)%num_proc] =  read_by_proc_offset[(rank-j+num_proc)%num_proc];
	}

	free(read_by_proc_offset);

	for(k=1; k<num_proc; k<<=1)
	{
		srank = (rank - k + num_proc) % num_proc;	//Rank to send to
		rrank = (rank + k) % num_proc;	//Rank to recv from
		int count = badCount(k, num_proc);
		recv_index = malloc(sizeof(int) * count);
		count = create_send_datatype_for_offsets(rank, num_proc, number_of_reads_by_procs,
				data_offsets, k, &dt_send, &recv_index);

		MPI_Type_commit(&dt_send);

		send_size_by_proc = malloc(count*sizeof(size_t));
		recv_size_by_proc = malloc(count*sizeof(size_t));

		send_total = get_send_size(rank, num_proc, number_of_reads_by_procs,
				&send_size_by_proc, count, k);

		MPI_Pack_size(1, dt_send, comm, &packsize);
		assert(packsize == (8 * send_total));

		MPI_Sendrecv(send_size_by_proc, count, MPI_LONG_LONG_INT, srank, 0,
				recv_size_by_proc, count, MPI_LONG_LONG_INT,
				rrank, 0, comm, MPI_STATUS_IGNORE);

		total = 0;

		for(m = 0; m < count; m++)
		{
			total += recv_size_by_proc[m];
		}
		size_t *interbuff_offset = calloc(total, sizeof(size_t));
		MPI_Sendrecv(MPI_BOTTOM, 1, dt_send, srank, 0,
				interbuff_offset, total, MPI_LONG_LONG_INT, rrank, 0, comm, &status);
		for ( m = 0; m < count; m++){
			// we free and allocate data_offsets
			// according to the recieve size
			if (data_offsets[recv_index[m]]){

				free(data_offsets[recv_index[m]]);
				data_offsets[recv_index[m]] = NULL;
				data_offsets[recv_index[m]] = malloc(sizeof(size_t)*(recv_size_by_proc[m]));
			}
		}
		size_t *tmp_var = interbuff_offset;

		for (m = 0; m < count; m++){

			memcpy(data_offsets[recv_index[m]], tmp_var, recv_size_by_proc[m] * sizeof(size_t));
			tmp_var += recv_size_by_proc[m];
			number_of_reads_by_procs[recv_index[m]] = recv_size_by_proc[m];

		}

		for (m = 0; m < count; m++){
			for (j = 0; j< number_of_reads_by_procs[recv_index[m]]; j++){
				//assert(data_offsets[recv_index[m]][j] != 0);
			}
		}
		MPI_Type_free(&dt_send);
		count = 0;

		// problem with the interbuff free !!!
		free(interbuff_offset);
		free(recv_index);
		free(recv_size_by_proc);
		free(send_size_by_proc);

	}
	free(data_offsets2);
}

void bruck_size(int rank, int num_proc, size_t local_readNum, size_t* number_of_reads_by_procs, int ** data_size, int *new_rank, int *new_size)
{
	MPI_Comm comm = COMM_WORLD;

	int k, m, j, srank, rrank;
	MPI_Datatype dt_send;
	size_t *recv_size_by_proc=NULL, *send_size_by_proc=NULL;
	int *recv_index=NULL;
	size_t total, send_total;
	int packsize;
	double time;

	int count;

	time = MPI_Wtime();

	for (m = 0; m < num_proc; m++){
		number_of_reads_by_procs[m] = 0;
	}

	for(m = 0; m < local_readNum; m++)
	{
		number_of_reads_by_procs[new_rank[m]]++;
	}

	int **data_size2 = malloc(sizeof(int *)*num_proc);
	//we initialize data_offsets
	for(m = 0; m < num_proc; m++){
		data_size[m] = NULL;
		data_size2[m] = calloc( number_of_reads_by_procs[m], sizeof(int));
	}


	size_t *read_by_proc = calloc(num_proc, sizeof(size_t));

	//we give values data_offsets
	for(m = 0; m < local_readNum; m++){

		// Phase one of bruck shift of the (rank+i)%size and (rank-i+size)%size
		data_size2[new_rank[m]][read_by_proc[new_rank[m]]] = new_size[m];
		read_by_proc[new_rank[m]]++;
	}

	for (j = 0; j < num_proc; j++){
		data_size[(rank+j)%num_proc] = data_size2[(rank-j+num_proc)%num_proc];
		number_of_reads_by_procs[(rank+j)%num_proc] =  read_by_proc[(rank-j+num_proc)%num_proc];
	}
	free(read_by_proc);

	for(k=1; k<num_proc; k<<=1)
	{
		srank = (rank - k + num_proc) % num_proc;	//Rank to send to
		rrank = (rank + k) % num_proc;				//Rank to recv from


		int count = badCount(k, num_proc);
		recv_index = malloc(sizeof(int) * count);

		count = create_send_datatype_for_size(rank, num_proc, number_of_reads_by_procs,
				data_size, k, &dt_send, &recv_index);

		MPI_Type_commit(&dt_send);

		send_size_by_proc = malloc(count * sizeof(size_t));
		recv_size_by_proc = malloc(count * sizeof(size_t));

		send_total = get_send_size(rank, num_proc, number_of_reads_by_procs,
				&send_size_by_proc, count, k);

		MPI_Pack_size(1, dt_send, comm, &packsize);

		MPI_Sendrecv(send_size_by_proc, count, MPI_LONG_LONG_INT, srank, 0,
				recv_size_by_proc, count, MPI_LONG_LONG_INT,
				rrank, 0, comm, MPI_STATUS_IGNORE);

		total = 0;

		for(m = 0; m < count; m++)
		{
			total += recv_size_by_proc[m];
		}

		int *interbuff_offset = malloc(total*sizeof(int));

		MPI_Sendrecv(MPI_BOTTOM, 1, dt_send, srank, 0,
				interbuff_offset, total, MPI_INT, rrank, 0, comm, MPI_STATUS_IGNORE);

		for ( m = 0; m < count; m++){
			// we free and allocate data_offsets
			// according to the recieve size
			if (data_size[recv_index[m]]){

				free(data_size[recv_index[m]]);
				data_size[recv_index[m]] = NULL;
				data_size[recv_index[m]] = malloc(sizeof(int)*(recv_size_by_proc[m]));
				}
		}

		int *tmp_var = interbuff_offset;

		for (m = 0; m < count; m++){

			memcpy(data_size[recv_index[m]], tmp_var, recv_size_by_proc[m] * sizeof(int));
			tmp_var += recv_size_by_proc[m];
			number_of_reads_by_procs[recv_index[m]] = recv_size_by_proc[m];

		}

		MPI_Type_free(&dt_send);

		count = 0;
		free(interbuff_offset);
		free(recv_index);
		free(recv_size_by_proc);
		free(send_size_by_proc);
	}
	free(data_size2);
}



void writeSam(
		int rank,
		char* output_dir,
		char* header,
		size_t local_readNum,
		size_t total_num_read,
		char *chrName,
		Read *chr,
		int total_num_proc,  //this is the number of proc in split communication
		MPI_Comm split_comm,
		int master_rank,
		char *file_name,
		MPI_File in,
		MPI_Info finfo,
		int compression_level,
		size_t *offset_dest_phase1,
		size_t *offset_source_phase1,
		int *read_size_phase1,
		int *dest_rank_phase1,
		// use when redistribute the reads according to original rank
		// when sorting of offset sources is done
		int *source_rank_phase1,
		char *data,
		size_t start_offset_in_file,
		size_t previous_local_readNum
		){


		/*
		 * This section is divided in 4 steps
		 *
		 * First step:
		 *
		 * To accelerate the writing the output are written in blocks. Each rank is going to write a block of contiguous reads.
		 * To do that we sort the offset destinations and gives a new rank to contigues read.
		 *
		 * Second step:
		 *
		 * Now read the reads. The reading is done in contiguous blocks.
		 * So we sort the offsets sources before reading the reads
		 *
		 * Third step:
		 *
		 * The shuffle of the read according to the offset destinations
		 * The shuffle is optimized with Bruck method.
		 *
		 * Fourth step:
		 *
		 * The writing.
		 *
		 */

		size_t j;
		size_t k;
		int ierr;

		//COMM_WORLD will become our new communicator
		COMM_WORLD = split_comm;

		int master_job_phase_1 = master_rank;
		int master_job_phase_2 = master_rank;

		MPI_Status status;

		/*
		 * phase 1 variables
		 */

		size_t* pbs_offset_dest_phase1;   // resized vector for all vectors in bitonic are same size
		size_t* pbs_offset_source_phase1;
		int* pbs_read_size_phase1;
		int* pbs_dest_rank_phase1;
		int* pbs_orig_rank_off_phase1;

		/*
		 * variables for first Bruck
		 */

		size_t *new_pbs_offset_source_phase1 	= NULL;
		size_t *new_pbs_offset_dest_phase1 		= NULL;
		int    *new_pbs_read_size_phase1 		= NULL;
		int    *new_pbs_dest_rank_phase1 		= NULL;
		int    *new_pbs_orig_rank_off_phase1 	= NULL;


		size_t *local_offset_destination_bruck;   		// = malloc(numItems*sizeof(size_t));
		size_t *local_offset_source_sorted_bruck;
		int *local_reads_sizes_sorted_bruck;
		int *local_reads_dest_rank_sorted_bruck;
		int *local_rank_source_offset_sorted_bruck;

		//variables for the writing part

		//the MPI datatype
		MPI_Datatype Datatype_Read_to_write;

		//variables for MPI writes and read
		MPI_File out;
		char* path;
		double time_count;
		//The data in which what we read will be kept
		char **data2 = malloc( total_num_proc * sizeof(char*));

		int dimensions = total_num_proc;

		size_t max_num_read = 0;
		MPI_Allreduce(&local_readNum, &max_num_read, 1, MPI_LONG_LONG_INT, MPI_MAX, split_comm);

		pbs_offset_dest_phase1   = calloc(max_num_read, sizeof(size_t));
		pbs_offset_source_phase1 = calloc(max_num_read, sizeof(size_t));
		pbs_orig_rank_off_phase1 = calloc(max_num_read, sizeof(size_t));
		pbs_read_size_phase1     = calloc(max_num_read, sizeof(int));
		pbs_dest_rank_phase1     = calloc(max_num_read, sizeof(int));

		size_t m=0;
		//we fill up the pbs vector
		for (m = 0; m < local_readNum; m++){
			pbs_offset_dest_phase1[m]   = offset_dest_phase1[m];
			pbs_offset_source_phase1[m] = offset_source_phase1[m];
			pbs_read_size_phase1[m]     = read_size_phase1[m];
			pbs_dest_rank_phase1[m]     = dest_rank_phase1[m];
			pbs_orig_rank_off_phase1[m] = source_rank_phase1[m];
		}

		free(offset_dest_phase1);
		free(offset_source_phase1);
		free(read_size_phase1);
		free(dest_rank_phase1);
		free(source_rank_phase1);

		/*
		 *
		 * CHECK OFFSET FOR DEBUG
		 *
		 * We check offset by reading the first character
		 * of the read in the SAM file
		 *
		for ( m = 0; m < local_readNum; m++){
			size_t offset_to_test = pbs_offset_source_phase1[m];
			char *buff_test = malloc(sizeof(char));
			buff_test[1] = 0;
			MPI_File_read_at(in, offset_to_test, buff_test, 1, MPI_CHAR, &status);
			assert( *buff_test == 'H' );
		}

		fprintf(stderr, "Rank %d :::::[WRITE][PHASE 1] check first char before bitonic 1 passed\n", rank, k);
		MPI_Barrier(COMM_WORLD);
		*/

		// we sort source offsets
		time_count = MPI_Wtime();
		ParallelBitonicSort2(
			COMM_WORLD,
			rank,
			dimensions,
			pbs_offset_source_phase1,
			pbs_read_size_phase1,
			pbs_dest_rank_phase1,
			pbs_offset_dest_phase1,
			pbs_orig_rank_off_phase1,
			max_num_read
		);

		MPI_Barrier(COMM_WORLD);

		if (rank == master_job_phase_1)
			fprintf(stderr, "Rank %d :::::[WRITE][BITONIC 2] Time spent sorting sources offsets = %f\n",
				rank,  MPI_Wtime() - time_count);

		/*
		 *
		 * FOR DEBUG
		 *
		 */

		for(j = 1; j < max_num_read; j++){
			assert(pbs_offset_source_phase1[j-1] <= pbs_offset_source_phase1[j]);
			assert(pbs_dest_rank_phase1[j] < dimensions);
			assert(pbs_orig_rank_off_phase1[j] < dimensions);
		}
		

		/*
		 *
		 * CHECK OFFSET
		 *
		 * We check offset by reading the first character
		 * of the read in the SAM file
		 *

		for ( m = 0; m < local_readNum; m++){

			if ( pbs_offset_source_phase1[m] > 0 ){
				size_t offset_to_test = pbs_offset_source_phase1[m];
				char *buff_test = malloc(sizeof(char));
				buff_test[1] = 0;
				MPI_File_read_at(in, offset_to_test, buff_test, 1, MPI_CHAR, &status);
				assert( *buff_test == 'H' );
			}

		}
		fprintf(stderr, "Rank %d :::::[WRITE][PHASE 1] check first char after bitonic 2 passed\n", rank, k);
		MPI_Barrier(COMM_WORLD);

		*/

		/*
		 * REMOVE ZERO PADDING BEFORE BRUCK
		 *
		 */


		size_t tmp2 = 0;
		for(j = 0; j < max_num_read; j++){ if (pbs_read_size_phase1[j] == 0) tmp2++;}

		MPI_Barrier(COMM_WORLD);

		size_t num_read_for_bruck = 0;
		local_readNum = max_num_read;

		if ( tmp2 < max_num_read ){
			//don't malloc!!
			new_pbs_offset_source_phase1 	= calloc( (local_readNum - tmp2), sizeof(size_t));
			new_pbs_offset_dest_phase1 		= calloc( (local_readNum - tmp2), sizeof(size_t));
			new_pbs_read_size_phase1 		= calloc( (local_readNum - tmp2), sizeof(int));
			new_pbs_dest_rank_phase1 		= calloc( (local_readNum - tmp2), sizeof(int));
			new_pbs_orig_rank_off_phase1 	= calloc( (local_readNum - tmp2), sizeof(int));

			for (j = 0; j < (local_readNum - tmp2); j++){
				new_pbs_offset_source_phase1[j] = pbs_offset_source_phase1[j+tmp2];
				new_pbs_read_size_phase1[j] 	= pbs_read_size_phase1[j+tmp2];
				new_pbs_dest_rank_phase1[j] 	= pbs_dest_rank_phase1[j+tmp2];
				new_pbs_orig_rank_off_phase1[j] = pbs_orig_rank_off_phase1[j+tmp2];
				new_pbs_offset_dest_phase1[j]	= pbs_offset_dest_phase1[j+tmp2];
			}

			/*
			 *
			 * FOR DEBUG
			 *
			for (j = 0; j < (local_readNum - tmp2); j++){
				assert(new_pbs_read_size_phase1[j]     != 0);
				assert(new_pbs_offset_source_phase1[j] != 0);
				assert(new_pbs_offset_dest_phase1[j]   != 0);
				assert(new_pbs_dest_rank_phase1[j]     < dimensions);
				assert(new_pbs_orig_rank_off_phase1[j] < dimensions);
			}
			*/
			num_read_for_bruck = local_readNum - tmp2;
		}

		if (tmp2 == max_num_read){

			size_t numItems = 0;

			new_pbs_offset_source_phase1    = malloc(numItems * sizeof(size_t));
			new_pbs_offset_dest_phase1      = malloc(numItems * sizeof(size_t));
			new_pbs_read_size_phase1        = malloc(numItems * sizeof(int));
			new_pbs_dest_rank_phase1        = malloc(numItems * sizeof(int));
			new_pbs_orig_rank_off_phase1    = malloc(numItems * sizeof(int));

			num_read_for_bruck = 0;
		}

		MPI_Barrier(COMM_WORLD);


		free(pbs_offset_source_phase1);
		free(pbs_read_size_phase1);
		free(pbs_dest_rank_phase1);
		free(pbs_orig_rank_off_phase1);
		free(pbs_offset_dest_phase1);
		/*
		 * Now we dipatch the vectors:
		 * 		read_size,
		 * 		dest_rank,
		 * 		offset_dest,
		 * 		offset_source
		 *
		 * according to their original ranks
		 * we do it with a Bruck
		 */

		int *new_local_reads_sizes_sorted_bruck 		= malloc(previous_local_readNum * sizeof(int));
		int *new_local_reads_dest_rank_sorted_bruck   	= malloc(previous_local_readNum * sizeof(int));
		size_t *new_local_offset_destination_bruck 		= malloc(previous_local_readNum * sizeof(size_t));
		size_t *new_local_offset_source_sorted_bruck	= malloc(previous_local_readNum * sizeof(size_t));

		int num_proc = dimensions;
		size_t *number_of_reads_by_procs = calloc( dimensions, sizeof(size_t));

		for(m = 0; m < num_read_for_bruck; m++){

			 number_of_reads_by_procs[new_pbs_orig_rank_off_phase1[m]]++;
		}

		size_t count6 = 0;
		for(m = 0; m < dimensions; m++){
				count6 += number_of_reads_by_procs[m];
		}
		assert( count6 == num_read_for_bruck );


		size_t **dest_offsets 			= malloc(sizeof(size_t *) * dimensions);
		size_t **local_source_offsets 	= malloc(sizeof(size_t *) * dimensions);

		int **read_size 				= malloc(sizeof(int *) * dimensions);
		int **dest_rank 				= malloc(sizeof(int *) * dimensions);


		/*
		 *
		 * CHECK OFFSET
		 *

		for ( m = 0; m < local_readNum; m++){

			if ( new_pbs_offset_source_phase1[m] > 0 ){
				size_t offset_to_test = new_pbs_offset_source_phase1[m];
				char *buff_test = malloc(sizeof(char));
				buff_test[1] = 0;
				MPI_File_read_at(in, offset_to_test, buff_test, 1, MPI_CHAR, &status);
				assert( *buff_test == 'H' );
			}
		}
		fprintf(stderr, "Rank %d :::::[WRITE][PHASE 1] check first char before bruck passed\n", rank, k);
		MPI_Barrier(COMM_WORLD);
		*/

		time_count = MPI_Wtime();
		bruckWrite2(
				rank,
				dimensions,
				count6,
				number_of_reads_by_procs,
				new_pbs_orig_rank_off_phase1,
				new_pbs_offset_dest_phase1,
				&dest_offsets,
				new_pbs_dest_rank_phase1,
				&dest_rank,
				new_pbs_offset_source_phase1,
				&local_source_offsets,
				new_pbs_read_size_phase1,
				&read_size
		);

		free(new_pbs_offset_dest_phase1);
		free(new_pbs_offset_source_phase1);
		free(new_pbs_read_size_phase1);
		free(new_pbs_dest_rank_phase1);
		free(new_pbs_orig_rank_off_phase1);

		/*
		 * Now get the offset source, destination and sizes
		 */
		MPI_Barrier(COMM_WORLD);
		if (rank == master_job_phase_1)
			fprintf(stderr, "Rank %d :::::[WRITE][BRUCK 2] Time spent in bruck  = %f\n", rank,  MPI_Wtime() - time_count);

		/*
		 * GET DATA AFTER BRUCK
		 *
		 */

		size_t count5=0;
		j=0;
		for(m = 0; m < num_proc; m++)
		{
			for(k = 0; k < number_of_reads_by_procs[m]; k++)
			{
				new_local_offset_source_sorted_bruck[k + j] 	= local_source_offsets[m][k];
				new_local_reads_dest_rank_sorted_bruck[k + j] 	= dest_rank[m][k];
				new_local_reads_sizes_sorted_bruck[k + j] 		= read_size[m][k];
				new_local_offset_destination_bruck[k + j] 		= dest_offsets[m][k];
			}
			free(local_source_offsets[m]);
			free(dest_rank[m]);
			free(read_size[m]);
			free(dest_offsets[m]);
			j += number_of_reads_by_procs[m];

		}
		assert( j == previous_local_readNum );

		for(k = 0; k < previous_local_readNum; k++){
			assert(new_local_offset_source_sorted_bruck[k] != 0);
		}


		time_count = MPI_Wtime();
		//init indices for qksort
		size_t *coord_index = malloc(previous_local_readNum * sizeof(size_t));

		for(j = 0; j < previous_local_readNum; j++){
			coord_index[j] = j;
		}

		base_arr2 = new_local_offset_source_sorted_bruck;
		qksort(coord_index, previous_local_readNum, sizeof(size_t), 0, previous_local_readNum - 1, compare_size_t);

		int *new_local_reads_sizes_sorted_bruck2 		= malloc(previous_local_readNum * sizeof(int));
		int *new_local_reads_dest_rank_sorted_bruck2   	= malloc(previous_local_readNum * sizeof(int));
		size_t *new_local_offset_destination_bruck2 	= malloc(previous_local_readNum * sizeof(size_t));
		size_t *new_local_offset_source_sorted_bruck2	= malloc(previous_local_readNum * sizeof(size_t));

		//We index data
		for(j = 0; j < previous_local_readNum; j++){

			new_local_reads_sizes_sorted_bruck2[j] 			= new_local_reads_sizes_sorted_bruck[coord_index[j]];
			new_local_reads_dest_rank_sorted_bruck2[j] 		= new_local_reads_dest_rank_sorted_bruck[coord_index[j]];
			new_local_offset_destination_bruck2 [j] 		= new_local_offset_destination_bruck[coord_index[j]];
			new_local_offset_source_sorted_bruck2[j]		= new_local_offset_source_sorted_bruck[coord_index[j]];
		}

		/*
		 *
		 * FOR DEBUG
		 *
		 *
		*/

		for(j = 0; j < previous_local_readNum - 1; j++){

			assert(new_local_offset_source_sorted_bruck2[j] < new_local_offset_source_sorted_bruck2[j+1]);
		}

		free(new_local_offset_source_sorted_bruck);
		free(new_local_reads_sizes_sorted_bruck);
		free(new_local_offset_destination_bruck);
		free(new_local_reads_dest_rank_sorted_bruck);
		free(coord_index);

		malloc_trim(0);

		if (rank == master_job_phase_2)
			fprintf(stderr, "Rank %d :::::[WRITE][LOCAL SORT] Time =  %f seconds\n", rank, MPI_Wtime() - time_count);

		/*
		 *
		 * FOR DEBUG
		 *
		 *
		*/
		for(k = 1; k < previous_local_readNum; k++)
		{
			assert((new_local_offset_source_sorted_bruck2[k] - start_offset_in_file) <= strlen(data));
			assert(new_local_offset_source_sorted_bruck2[k-1]    <= new_local_offset_source_sorted_bruck2[k]);
			assert(new_local_offset_source_sorted_bruck2[k] 	 != 0);
			assert(new_local_reads_dest_rank_sorted_bruck2[k] 	 < dimensions);
			assert(new_local_offset_destination_bruck2[k] 		 != 0);
			assert(new_local_reads_sizes_sorted_bruck2[k] 		 != 0);
		}
		 

		if (dest_rank != NULL)
			free(dest_rank);
		if (read_size != NULL)
			free(read_size);
		if (dest_offsets != NULL)
			free(dest_offsets);
		if (local_source_offsets != NULL)
			free(local_source_offsets);

		/*
		 *
		 * PACKING OF THE DATA
		 *
		 */

		size_t i;
		size_t new_data_sz = 0;
		char *data_pack;
		char *tmp_tab;

		//we compute the size of data_pack
		for (k = 0; k < previous_local_readNum; k++){
			new_data_sz += new_local_reads_sizes_sorted_bruck2[k];
		}

		data_pack = malloc(new_data_sz +1);
		data_pack[new_data_sz] = 0;

		char *q = data;
		char *p = data_pack;
		size_t offset_in_data = 0;
		int pos = 0;

		MPI_Barrier(COMM_WORLD);

		//we copy elements from data in data_pack
		for (k=0; k < (previous_local_readNum); k++){
			pos = 0;
			offset_in_data = new_local_offset_source_sorted_bruck2[k] - start_offset_in_file;
			q = data + offset_in_data;
			while (*q && (pos < new_local_reads_sizes_sorted_bruck2[k])) { *p=*q; q++; p++; pos++; }
		}

		free(new_local_offset_source_sorted_bruck2);
		int res;

		/*
		 * We unpack in a loop the same way
		 */

		MPI_Datatype dt_data;
		time_count = MPI_Wtime();
		//The data in which what we read will be kept
		// we compute the size of
		// data for each rank and we put it in
		// buffs and buffs_by_proc is
		// the size of the buffer to send

		size_t *buffs_by_procs2 = calloc( dimensions, sizeof(size_t));
		size_t *buffs_by_procs  = calloc( dimensions, sizeof(size_t));

		for(m = 0; m < num_proc; m++) number_of_reads_by_procs[m] = 0;

		for(m = 0; m < previous_local_readNum; m++)
		{
			buffs_by_procs2[new_local_reads_dest_rank_sorted_bruck2[m]] += new_local_reads_sizes_sorted_bruck2[m];
			number_of_reads_by_procs[new_local_reads_dest_rank_sorted_bruck2[m]]++;
		}

		for(m = 0; m < num_proc; m++)
		{
			buffs_by_procs[(rank + m)%num_proc] = buffs_by_procs2[(rank - m + num_proc)%num_proc];
		}

		free(buffs_by_procs2);

		//Allocate data and initialization
		for(m = 0; m < num_proc; m++)
		{
			data2[m] = (char*)malloc(buffs_by_procs[m]*sizeof(char) + 1);
			data2[m][buffs_by_procs[m]] = 0;
		}

		//Variable for datatype struct
		MPI_Aint *indices 		= malloc(previous_local_readNum * sizeof(MPI_Aint));
		int *blocklens    		= malloc(previous_local_readNum * sizeof(int));
		MPI_Datatype *oldtypes 	= malloc(previous_local_readNum * sizeof(MPI_Datatype));

		MPI_Aint adress_to_write_in_data_by_element[num_proc];
		for(i = 0; i < num_proc; i++){
			MPI_Get_address(data2[(rank-i+num_proc)%num_proc], &adress_to_write_in_data_by_element[(rank+i)%num_proc]);
		}

		for(i = 0; i < previous_local_readNum; i++){
			indices[i] = adress_to_write_in_data_by_element[new_local_reads_dest_rank_sorted_bruck2[i]];
			assert (indices[i] != (MPI_Aint)NULL);
			adress_to_write_in_data_by_element[new_local_reads_dest_rank_sorted_bruck2[i]] += new_local_reads_sizes_sorted_bruck2[i];
			blocklens[i] = new_local_reads_sizes_sorted_bruck2[i];
			oldtypes[i] = MPI_CHAR;
		 }

		//Create struct
		MPI_Type_create_struct(previous_local_readNum, blocklens, indices, oldtypes, &dt_data);
		MPI_Type_commit(&dt_data);
		pos=0;
		res = MPI_Unpack(data_pack, new_data_sz, &pos, MPI_BOTTOM, 1, dt_data, COMM_WORLD);
		assert(res == MPI_SUCCESS);
		MPI_Type_free(&dt_data);

		free(data_pack);
		free(blocklens);
		free(indices);
		free(oldtypes);

		malloc_trim(0);
		if (rank == master_job_phase_2)
			fprintf(stderr, "Rank %d :::::[WRITE][DATA PACK] Time =  %f seconds\n", rank, MPI_Wtime() - time_count);

		/*
		 * 	In this part we are going to send the data
		 * 	buffer according to the rank it belong
		 * 	and in the sorted order
		 *
		 *  variable we have are:
		 *
		 *  	1) data the buffer = hold the reads
		 *  	2) new_rank_sorted_phase2 vector = hold the rank of the reads
		 *		3) new_offset_dest_sorted_phase2 = hold the offsets of the in the destination file
		 *		4) new_read_size_sorted_phase2 = hold the size of the reads
		 *
		 *
		 *
		 *	The strategy is :
		 *		1) we loop the rank
		 *		2) each rank send to the target rank how much data it's going to send. In order to prepare buffers
		 *		3) for each rank we create datatype of buffered data, read size, and offset
		 *		4) the taget rank recieve the buffered data, the reads size, and the offset vector
		 */

		/****************************
		 * 	BEGIN BRUCK PHASE     	*
		 *****************************/

		size_t **data_offsets = malloc(sizeof(size_t *) * num_proc);
		int **data_size 	  = malloc(sizeof(int *) * num_proc);

		time_count = MPI_Wtime();

		count6 = 0;
		for(m = 0; m < dimensions; m++){
				count6 += number_of_reads_by_procs[m];
		}
		assert( count6 == previous_local_readNum );

		bruckWrite(
				rank,
				num_proc,
				previous_local_readNum,
				number_of_reads_by_procs,
				new_local_reads_dest_rank_sorted_bruck2,
				buffs_by_procs,
				&data2,
				new_local_offset_destination_bruck2,
				&data_offsets,
				new_local_reads_sizes_sorted_bruck2,
				&data_size
		);

		if (rank == master_job_phase_1)
			fprintf(stderr, "Rank %d :::::[WRITE][BRUCK] Time   = %f s \n", rank, MPI_Wtime() - time_count);

		for(m = 0; m < num_proc; m++)
		{
			size_t buff_test  = strlen(data2[m]);
			size_t buff_test2 = 0;
			size_t u = 0;
			for (u = 0; u < number_of_reads_by_procs[m]; u++)
				buff_test2 += data_size[m][u];

			assert (buff_test2 == buff_test);
		}

		free(buffs_by_procs);
		free(new_local_reads_sizes_sorted_bruck2);
		free(new_local_offset_destination_bruck2);
		free(new_local_reads_dest_rank_sorted_bruck2);
		malloc_trim(0);

		size_t sum_num_of_reads = 0;
		for (m = 0; m < num_proc; m++){
			sum_num_of_reads += number_of_reads_by_procs[m];
		}
		assert(previous_local_readNum == sum_num_of_reads);

		size_t *new_offset_dest_index_phase3	= malloc(sizeof(size_t) * previous_local_readNum);
		char  **data_reads_to_sort 				= malloc(previous_local_readNum * sizeof(char*));
		int *data_size_to_sort 					= malloc(previous_local_readNum * sizeof(int));
		size_t *data_offsets_to_sort 			= malloc(previous_local_readNum * sizeof(size_t));

		/*
		 * GET DATA AFTER BRUCK
		 *
		 */

		j=0;
		for(m = 0; m < num_proc; m++)
		{
			int i = 0;
			for(k = 0; k < number_of_reads_by_procs[m]; k++)
			{
				data_reads_to_sort[k + j] = &(data2[m][i]);
				i += data_size[m][k];

				data_size_to_sort[k + j] = data_size[m][k];
				data_offsets_to_sort[k + j] = data_offsets[m][k];
			}
			free(data_size[m]);
			free(data_offsets[m]);
			j += number_of_reads_by_procs[m];
		}
		assert(j == previous_local_readNum);

		/*
		 *
		 * FOR DEBUG
		 *
		 *

		for ( k = 0; k < previous_local_readNum; k++){
			assert (data_size_to_sort != 0);
			assert( data_offsets_to_sort != 0);
		}
		 */

		if (data_offsets != NULL)
			free(data_offsets);
		if (data_size != NULL)
			free(data_size);

		free(number_of_reads_by_procs);

		/*
		 * SORT LOCALY OFFSET DESTINATION BEFORE WRITING
		 *
		 */

		for (k = 0; k < previous_local_readNum; k++){
			new_offset_dest_index_phase3[k] = k;
		}
		//task SORT OFFSET DESTINATION
		//now we sort new_offset_dest_phase2


		base_arr2 = data_offsets_to_sort;
		qksort(new_offset_dest_index_phase3, previous_local_readNum, sizeof(size_t), 0, previous_local_readNum - 1, compare_size_t);


		size_t* offsets_sorted = malloc(sizeof(size_t) * previous_local_readNum);
		for (k = 0; k < previous_local_readNum; k++){
			offsets_sorted[k] = data_offsets_to_sort[new_offset_dest_index_phase3[k]];
		}


		for (k = 0; k < previous_local_readNum - 1; k++){
			assert (offsets_sorted[k] < offsets_sorted[k+1]);
		}


		free(data_offsets_to_sort);

		size_t size_t_buffer_uncompressed = 0;
		for(k = 0; k < previous_local_readNum; k++){
			size_t_buffer_uncompressed += data_size_to_sort[new_offset_dest_index_phase3[k]];
		}

		char *char_buff_uncompressed = malloc(size_t_buffer_uncompressed * sizeof(char) + 1);
		char_buff_uncompressed[size_t_buffer_uncompressed] = 0;
		char *p1 = char_buff_uncompressed;
		size_t q1 = 0;

		for(k = 0; k < previous_local_readNum; k++){
			while ( q1 < data_size_to_sort[new_offset_dest_index_phase3[k]]){
				*p1++ = *data_reads_to_sort[new_offset_dest_index_phase3[k]]++;
				q1++;
			}
			q1=0;
		}

		free(new_offset_dest_index_phase3);
		free(data_reads_to_sort);

		//fprintf (stderr, "[WRITE] rank = %d :: char_buff_uncompressed = %s \n", rank, char_buff_uncompressed);

		/*
		 * COMPRESSION PART
		 *
		 */

		time_count = MPI_Wtime();

		BGZF *fp;
		fp = calloc(1, sizeof(BGZF));
		int block_length = MAX_BLOCK_SIZE;
		int bytes_written;
		int length = strlen(char_buff_uncompressed);

		fp->open_mode = 'w';
		fp->uncompressed_block_size = MAX_BLOCK_SIZE;
		fp->uncompressed_block = malloc(MAX_BLOCK_SIZE);
		fp->compressed_block_size = MAX_BLOCK_SIZE;
		fp->compressed_block = malloc(MAX_BLOCK_SIZE);
		fp->cache_size = 0;
		fp->cache = kh_init(cache);
		fp->block_address = 0;
		fp->block_offset = 0;
		fp->block_length = 0;
		fp->compress_level = compression_level < 0? Z_DEFAULT_COMPRESSION : compression_level; // Z_DEFAULT_COMPRESSION==-1

		if (fp->compress_level > 9) fp->compress_level = Z_DEFAULT_COMPRESSION;

		const bgzf_byte_t *input = (void *)char_buff_uncompressed;
		int compressed_size = 0;

		if (fp->uncompressed_block == NULL)
		   fp->uncompressed_block = malloc(fp->uncompressed_block_size);

		input = (void *)char_buff_uncompressed;
		block_length = fp->uncompressed_block_size;
		bytes_written = 0;
		uint8_t *compressed_buff =  malloc(strlen(char_buff_uncompressed) * sizeof(uint8_t));

		if (rank == master_job_phase_2)
			fprintf(stderr, "rank %d :::: start loop compression \n", rank);

		while (bytes_written < length) {
			int copy_length = bgzf_min(block_length - fp->block_offset, length - bytes_written);
			bgzf_byte_t* buffer = fp->uncompressed_block;
			memcpy(buffer + fp->block_offset, input, copy_length);
			fp->block_offset += copy_length;
			input += copy_length;
			bytes_written += copy_length;
			//if (fp->block_offset == block_length) {
				//we copy in a temp buffer
				while (fp->block_offset > 0) {
					int block_length;
					block_length = deflate_block(fp, fp->block_offset);

					//is it necessary?
					//if (block_length < 0) break;

					// count = fwrite(fp->compressed_block, 1, block_length, fp->file);
					// we replace the fwrite with a memcopy
					memcpy(compressed_buff + compressed_size, fp->compressed_block, block_length);
					compressed_size +=block_length;
					fp->block_address += block_length;
				}
			//}
		}

		BGZF *fp_header;
		fp_header = calloc(1, sizeof(BGZF));
		uint8_t *compressed_header = NULL;
		int compressed_size_header = 0;

		if (rank == 0) {

			int block_length = MAX_BLOCK_SIZE;
			int bytes_written;
			int length = strlen(header);

			fp_header->open_mode = 'w';
			fp_header->uncompressed_block_size = MAX_BLOCK_SIZE;
			fp_header->uncompressed_block = malloc(MAX_BLOCK_SIZE);
			fp_header->compressed_block_size = MAX_BLOCK_SIZE;
			fp_header->compressed_block = malloc(MAX_BLOCK_SIZE);
			fp_header->cache_size = 0;
			fp_header->block_address = 0;
			fp_header->block_offset = 0;
			fp_header->block_length = 0;
			fp_header->compress_level = compression_level < 0? Z_DEFAULT_COMPRESSION : compression_level; // Z_DEFAULT_COMPRESSION==-1

			if (fp_header->compress_level > 9) fp_header->compress_level = Z_DEFAULT_COMPRESSION;


			const bgzf_byte_t *input = (void *)header;


			if (fp_header->uncompressed_block == NULL)
				fp_header->uncompressed_block = malloc(fp_header->uncompressed_block_size);

			input = (void *)header;
			block_length = fp_header->uncompressed_block_size;
			bytes_written = 0;
			compressed_header =  malloc(strlen(char_buff_uncompressed) * sizeof(uint8_t));

			fprintf(stderr, "rank %d :::: start loop compression \n", rank);

			while (bytes_written < length) {
				int copy_length = bgzf_min(block_length - fp_header->block_offset, length - bytes_written);
				bgzf_byte_t* buffer = fp_header->uncompressed_block;
				memcpy(buffer + fp_header->block_offset, input, copy_length);
				fp_header->block_offset += copy_length;
				input += copy_length;
				bytes_written += copy_length;
				//if (fp->block_offset == block_length) {
				//we copy in a temp buffer
				while (fp_header->block_offset > 0) {
					int block_length;
					block_length = deflate_block(fp_header, fp_header->block_offset);

					//is it necessary?
					//if (block_length < 0) break;

					// count = fwrite(fp->compressed_block, 1, block_length, fp->file);
					// we replace the fwrite with a memcopy
					memcpy(compressed_header + compressed_size_header, fp_header->compressed_block, block_length);
					compressed_size_header +=block_length;
					fp_header->block_address += block_length;
				}
				//}
			}
		}

	        kh_destroy(cache, fp->cache);
		free(char_buff_uncompressed);
		size_t compSize = compressed_size;

		if (rank == master_job_phase_2)
				fprintf(stderr, "Rank %d :::::[WRITE][COMPRESSION] Time for compressing %f seconds \n",
					rank, MPI_Wtime()-time_count);



		/*
		 * We write results of compression
		 */

		size_t write_offset = 0;

		MPI_Offset * y  = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
		MPI_Offset * y2 = (MPI_Offset *) calloc(num_proc+1, sizeof(MPI_Offset));

		MPI_Gather(&compSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);

		//now we make a cumulative sum
		int i1 = 0;

		if (rank ==0){
			for (i1 = 1; i1 < (num_proc + 1); i1++) {
				y2[i1] = y[i1-1];
				}

			for (i1 = 1; i1 < (num_proc +1); i1++) {
				y2[i1] = y2[i1-1] + y2[i1];
			}

			for (i1 = 0; i1 < (num_proc +1); i1++) {
				y2[i1] = y2[i1] + write_offset + compressed_size_header;
			}

		}

		MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);
		// we create the path where to write for collective write
		path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
		sprintf(path, "%s/%s.gz", output_dir, chrName);


		//task FINE TUNING FINFO FOR WRITING OPERATIONS

		MPI_Info_set(finfo,"striping_factor","12");
		//MPI_Info_set(finfo,"striping_unit","1610612736"); //1G striping
		MPI_Info_set(finfo,"striping_unit","268435456"); //256 Mo

		MPI_Info_set(finfo,"nb_proc","12");
		MPI_Info_set(finfo,"cb_nodes","12");
		MPI_Info_set(finfo,"cb_block_size","4194304"); /* 4194304 = 4 MBytes - should match FS block size */
		//MPI_Info_set(finfo,"cb_buffer_size","1610612736"); /* 128 MBytes (Optional) */


		ierr = MPI_File_open(COMM_WORLD, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);

		if (ierr) {
			fprintf(stderr, "Rank %d failed to open %s.\nAborting.\n\n", rank, path);
			MPI_Abort(COMM_WORLD, ierr);
			exit(2);
		}
		else{
			if(!rank)fprintf(stderr, "Rank %d :::::[WRITE] %s.bam successfully opened\n", rank, chrName);
		}

		time_count = MPI_Wtime();

		if (rank == master_job_phase_2 ) {
			fprintf(stderr, "Proc rank %d ::: we write the header \n", rank);
			MPI_File_write(out, compressed_header, compressed_size_header, MPI_BYTE, MPI_STATUS_IGNORE);
		}
		free(compressed_header);

		MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
		MPI_File_write_all(out, compressed_buff, (size_t)compSize, MPI_BYTE, &status);

		//task FINE TUNING FINFO BACK TO READING OPERATIONS
		//MPI_Info_set(finfo,"striping_factor","12");
		//MPI_Info_set(finfo,"striping_unit","2684354560"); //1G striping

		//MPI_Info_set(finfo,"nb_proc","128");
		//MPI_Info_set(finfo,"cb_nodes","128");
		//MPI_Info_set(finfo,"cb_block_size","2684354560"); /* 4194304 MBytes - should match FS block size */
		//MPI_Info_set(finfo,"cb_buffer_size","2684354560"); /* 128 MBytes (Optional) */

		//free(buff_compressed);
		free(compressed_buff);

		if (rank == master_job_phase_2)
			fprintf(stderr, "Rank %d :::::[WRITE][WRITING BGZF] Time for chromosome %s writing %f seconds\n",
					rank, chrName, MPI_Wtime()-time_count);

		free(fp->uncompressed_block);
		free(fp->compressed_block);
		free_cache(fp);
		free(fp);

		if (rank == 0){
			free(fp_header->uncompressed_block);
			free(fp_header->compressed_block);
			free_cache(fp_header);
		}
		free(fp_header);


		MPI_File_close(&out);
		free(path);
		for(m = 0; m < num_proc; m++)
		{
			if (data2[m]) free(data2[m]);
		}
		if (data2) free(data2);
		free(offsets_sorted);
		free(data_size_to_sort);
		free(y);
		free(y2);
		malloc_trim(0);
}

void writeSam_discordant_and_unmapped(
		int split_rank,
		char* output_dir,
		char* header,
		size_t local_readNum,
		char* chrName,
		Read* chr,
		int split_size,
		MPI_Comm split_comm,
		char *file_name,
		MPI_File in,
		MPI_Info finfo,
		int compression_level,
		char* data,
		size_t start_offset_in_file,
		int write_sam
		){

	/*
	 * task: writeSam_unmapped write unmapped reads
	 *
	 */

	MPI_Status status;
	size_t j;
	size_t k;
	int ierr;
	size_t dataSize;

	// vector use in the reading and writing part
	size_t *offset_source_index = (size_t*)malloc(local_readNum*sizeof(size_t));
	size_t *offset_source_sorted = (size_t*)malloc(local_readNum*sizeof(size_t));
	int *read_size_unsorted = (int *) malloc(local_readNum * sizeof(int));
	int *read_size_sorted = (int *)malloc(local_readNum * sizeof(int));
	offset_source_index[0] = 0;

	//variables for MPI writes and read
	//MPI_File in, out;
	MPI_File out;
	char* path;
	double start, finish, io_time, time_count;

	//size_t offset_source[local_readNum];
	size_t *offset_source_unsorted = (size_t*)malloc(local_readNum*sizeof(size_t));
	offset_source_unsorted[0] = 0;

	//int size_source[local_readNum];

	read_size_unsorted[0] = 0;
	dataSize = 0;

	int master_job = 0;
	double start_phase2, finish_phase2;

	//we initialize offset source and size_source
	for(j = 0; j < local_readNum; j++){
		 read_size_unsorted[j] = 0;
		 offset_source_unsorted[j] = 0;
	}

	char *char_buff_uncompressed = malloc(1024*sizeof(char));
	assert(char_buff_uncompressed != 0);

	size_t *offset_in_data = malloc(local_readNum*sizeof(size_t));
	assert(offset_in_data != 0);

	uint8_t *char_buff_compressed = malloc(1024*sizeof(uint8_t));
	assert(char_buff_compressed != 0);
	char_buff_compressed[0] = 0;

	uint8_t *compressed_header = malloc(1024*sizeof(uint8_t));
	assert(compressed_header != 0);
	compressed_header[0] = 0;

	//COMM_WORLD will become our new communicator
	COMM_WORLD = split_comm;

	//variables for the writing part
	//the MPI datatype
	MPI_Datatype Datatype_Read_to_write;
	//variables for MPI writes and read

	//first we parse the chr structure
	//and get information of the offset in the sourcefile
	//offsets are then translated into the data buffer offset
	for(j = 0; j < local_readNum; j++){
		//offset is the read size
		offset_source_unsorted[j] = chr->offset_source_file;
		read_size_unsorted[j] = (int)chr->offset; //read length
		dataSize += chr->offset;
		chr = chr->next;
	}

	for(j = 1; j < local_readNum; j++){
		assert(offset_source_unsorted[j-1] < offset_source_unsorted[j]);
	}

	finish_phase2 = MPI_Wtime();
	io_time = finish_phase2 - start_phase2;
	//step 1 :: we sort the vector of input offset
	//we sort the offset_sources
	start = MPI_Wtime();
	for(j = 0; j < local_readNum; j++){
		offset_source_index[j] = j;
	}


	//previous version of local sort with output of permutation
	base_arr2 = offset_source_unsorted;
	//new version of the local sort
	qksort(offset_source_index, local_readNum, sizeof(size_t), 0, local_readNum - 1, compare_size_t);

	/*
	 * reorder the offset and the size
	 * of the reads according to sort
	 */

	for(j = 0; j < local_readNum; j++){
		offset_source_sorted[j] = offset_source_unsorted[offset_source_index[j]];
		read_size_sorted[j] = read_size_unsorted[offset_source_index[j]];
	}

	for(j = 1; j < local_readNum; j++){
		assert(offset_source_sorted[j-1] <= offset_source_sorted[j]);
		assert(read_size_sorted[j] !=0 );
	}

	free(offset_source_unsorted);
	free(read_size_unsorted);
	MPI_Barrier(split_comm);

	int m;
	size_t i;
	size_t new_data_sz = 0;
	//we compute the size of data_pack
	//and offset in data

	for (k = 0; k < local_readNum; k++){
		new_data_sz += read_size_sorted[k];
		offset_in_data[k] = offset_source_sorted[k] - start_offset_in_file;
	}

	char_buff_uncompressed = realloc(char_buff_uncompressed, new_data_sz +1);
	char_buff_uncompressed[new_data_sz] = 0;

	char *q = data;
	char *p = char_buff_uncompressed;
	int pos = 0;
	//we compute the new offset of reads in data buffer
	//we remove the start offset in the file

	size_t total_copy=0;
	//we copy elements from data in data_pack
	for (k=0; k < local_readNum; k++){
		pos = 0;
		q = data + offset_in_data[k];
		while (*q && (pos < read_size_sorted[k])) {*p=*q;q++;p++;pos++;}
	}

	BGZF *fp;
	fp = calloc(1, sizeof(BGZF));
	int compress_level = compression_level;
	int block_length = MAX_BLOCK_SIZE;
	int bytes_written;
	int length = strlen(char_buff_uncompressed);

	fp->open_mode = 'w';
	fp->uncompressed_block_size = MAX_BLOCK_SIZE;
	fp->uncompressed_block = malloc(MAX_BLOCK_SIZE);
	fp->compressed_block_size = MAX_BLOCK_SIZE;
	fp->compressed_block = malloc(MAX_BLOCK_SIZE);
	fp->cache_size = 0;
	fp->cache = kh_init(cache);
	fp->block_address = 0;
	fp->block_offset = 0;
	fp->block_length = 0;
	fp->compress_level = compress_level < 0? Z_DEFAULT_COMPRESSION : compress_level; // Z_DEFAULT_COMPRESSION==-1

	if (fp->compress_level > 9) fp->compress_level = Z_DEFAULT_COMPRESSION;

	const bgzf_byte_t *input = (void *)char_buff_uncompressed;
	size_t compressed_size = 0;

	if (fp->uncompressed_block == NULL)
	   fp->uncompressed_block = malloc(fp->uncompressed_block_size);

	input = (void *)char_buff_uncompressed;
	block_length = fp->uncompressed_block_size;
	bytes_written = 0;
	char_buff_compressed =  realloc(char_buff_compressed, (strlen(char_buff_uncompressed)+1) * sizeof(uint8_t));
	assert(char_buff_compressed != 0);
	char_buff_compressed[strlen(char_buff_uncompressed)]=0;

	time_count = MPI_Wtime();
	while (bytes_written < length) {
		int copy_length = bgzf_min(block_length - fp->block_offset, length - bytes_written);
		bgzf_byte_t* buffer = fp->uncompressed_block;
		memcpy(buffer + fp->block_offset, input, copy_length);
		fp->block_offset += copy_length;
		input += copy_length;
		bytes_written += copy_length;
		//if (fp->block_offset == block_length) {
			//we copy in a temp buffer
			while (fp->block_offset > 0) {
				int block_length;
				block_length = deflate_block(fp, fp->block_offset);
					//is it necessary?
				//if (block_length < 0) break;
					// count = fwrite(fp->compressed_block, 1, block_length, fp->file);
				// we replace the fwrite with a memcopy
				memcpy(char_buff_compressed + compressed_size, fp->compressed_block, block_length);
				compressed_size +=block_length;
				fp->block_address += block_length;
			}
		//}
	} //end data compression

	//we compress the neader
	BGZF *fp_header;
	fp_header = calloc(1, sizeof(BGZF));
	//uint8_t *compressed_header = NULL;
	size_t compressed_size_header = 0;

	if (split_rank == 0) {

		int compress_level = compression_level;
		int block_length = MAX_BLOCK_SIZE;
		int bytes_written;
		int length = strlen(header);

		fp_header->open_mode = 'w';
		fp_header->uncompressed_block_size = MAX_BLOCK_SIZE;
		fp_header->uncompressed_block = malloc(MAX_BLOCK_SIZE);
		fp_header->compressed_block_size = MAX_BLOCK_SIZE;
		fp_header->compressed_block = malloc(MAX_BLOCK_SIZE);
		fp_header->cache_size = 0;
		fp_header->block_address = 0;
		fp_header->block_offset = 0;
		fp_header->block_length = 0;
		fp_header->compress_level = compress_level < 0? Z_DEFAULT_COMPRESSION : compress_level; // Z_DEFAULT_COMPRESSION==-1

		if (fp_header->compress_level > 9) fp_header->compress_level = Z_DEFAULT_COMPRESSION;

		const bgzf_byte_t *input = (void *)header;

		if (fp_header->uncompressed_block == NULL)
			fp_header->uncompressed_block = malloc(fp_header->uncompressed_block_size);

		input = (void *)header;
		block_length = fp_header->uncompressed_block_size;
		bytes_written = 0;
		compressed_header =  realloc(compressed_header, (strlen(char_buff_uncompressed) + 1) * sizeof(uint8_t));
		compressed_header[strlen(char_buff_uncompressed)] = 0;

		while (bytes_written < length) {

			int copy_length = bgzf_min(block_length - fp_header->block_offset, length - bytes_written);
			bgzf_byte_t* buffer = fp_header->uncompressed_block;
			memcpy(buffer + fp_header->block_offset, input, copy_length);
			fp_header->block_offset += copy_length;
			input += copy_length;
			bytes_written += copy_length;
			//if (fp->block_offset == block_length) {
			//we copy in a temp buffer
			while (fp_header->block_offset > 0) {
				int block_length;
				block_length = deflate_block(fp_header, fp_header->block_offset);
				//is it necessary?
				//if (block_length < 0) break;
				// count = fwrite(fp->compressed_block, 1, block_length, fp->file);
				// we replace the fwrite with a memcopy
				memcpy(compressed_header + compressed_size_header, fp_header->compressed_block, block_length);
				compressed_size_header +=block_length;
				fp_header->block_address += block_length;
			}
			//}
		}
	} //end header compression

        kh_destroy(cache, fp->cache);

	MPI_Barrier(split_comm);
	size_t write_offset = 0;

	MPI_Offset * y = (MPI_Offset *) calloc(split_size, sizeof(MPI_Offset));
	MPI_Offset * y2 = (MPI_Offset *) calloc(split_size+1, sizeof(MPI_Offset));

	MPI_Gather(&compressed_size, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, split_comm);

	//we make a cumulative sum
	int i1 = 0;

	if (split_rank ==0){
		for (i1 = 1; i1 < (split_size + 1); i1++) {
			y2[i1] = y[i1-1];
		}

		for (i1 = 1; i1 < (split_size +1); i1++) {
			y2[i1] = y2[i1-1] + y2[i1];
		}

		for (i1 = 0; i1 < (split_size +1); i1++) {
			y2[i1] = y2[i1] + write_offset + compressed_size_header;
		}
	}

	//do a gather in replacement of the the ring pass
	MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, split_comm);


	// we create the path where to write for collective write
	path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
	sprintf(path, "%s/%s.gz", output_dir, chrName);

	ierr = MPI_File_open(split_comm, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);

	if (ierr) {
		fprintf(stderr, "Rank %d :::::[WRITE] failed to open %s.\nAborting.\n\n", split_rank, path);
		MPI_Abort(split_comm, ierr);
		exit(2);
	}

	time_count = MPI_Wtime();

	if (split_rank == 0 ) {
		MPI_File_write(out, compressed_header, compressed_size_header, MPI_BYTE, MPI_STATUS_IGNORE);
		//we update write _header
	}

	MPI_Barrier(split_comm);
	//task WRITING OPERATIONS FOR UNMAPPED READS
	MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
	MPI_File_write(out, char_buff_compressed, (size_t)compressed_size, MPI_BYTE, &status);
	if (split_rank == master_job)
		fprintf(stderr, "Rank %d :::::[WRITE] Time for chromosome %s writing %f seconds\n", split_rank, chrName, MPI_Wtime()-time_count);

	free(compressed_header);
	free(fp->uncompressed_block);
	free(fp->compressed_block);
	free_cache(fp);
	free(fp);

	if (split_rank == 0){
		free(fp_header->uncompressed_block);
		free(fp_header->compressed_block);
		free_cache(fp_header);
	}

	free(fp_header);
	free(offset_in_data);
	MPI_File_close(&out);
	free(char_buff_compressed);
	free(char_buff_uncompressed);
	free(read_size_sorted);
	free(offset_source_index);
	free(offset_source_sorted);
	free(y);
	free(y2);
	free(path);
	malloc_trim(0);
}


void writeSam_any_dim(
		int dimensions,
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
		char *data,
		size_t start_offset_in_file){


	/*
	 *
	 * This section is divided in 4 steps
	 *
	 * First step:
	 *
	 * To accelerate the writing the output are written in blocks. Each rank is going to write a block of contiguous reads.
	 * To do that we sort the offset destinations and gives a new rank to contigues read.
	 *
	 * Second step:
	 *
	 * We extract the reads from the local buffer in memory.
	 * To do so we sort the offsets sources before reading the reads.
	 *
	 * Third step:
	 *
	 * The shuffle of the read according to the rank destinations
	 * The shuffle is optimized with Bruck method.
	 *
	 * Fourth step:
	 *
	 * After dispatching the reads. The output offsets are sorted
	 * and the reads are written in contigous blocks.
	 *
	 */

	size_t j;
	size_t k;
	int ierr;

	//COMM_WORLD will become our new communicator
	COMM_WORLD = split_comm;

	int master_job_phase_1 = master_rank;
	int master_job_phase_2 = master_rank;

	MPI_Status status;

	int *all_rank_to_send					= NULL;
	int *all_read_size_to_send				= NULL;
	size_t *all_offset_dest_file_to_send	= NULL;
	size_t *all_offset_source_file_to_send	= NULL;

	/*
	 * phase 1 variables
	 */

	int *all_read_size_phase1_to_send				= NULL;
	int *all_rank_phase1_to_send					= NULL;
	int *new_read_size_phase1 						= malloc(local_readNum* sizeof(int));
	int *new_rank_phase1 							= malloc(local_readNum* sizeof(int));

	size_t* all_offset_source_file_phase1_to_send	= NULL;
	size_t *new_offset_dest_phase1 					= malloc(local_readNum* sizeof(size_t));
	size_t *new_offset_source_phase1 				= malloc(local_readNum* sizeof(size_t));
	size_t *pbs_local_dest_offset					= NULL;
	size_t *pbs_local_dest_offset_index				= NULL;
	size_t *all_offset_dest_sorted_index_phase1		= NULL;
	size_t *all_offset_dest_file_to_send_phase1	    = NULL;

	/*
	 * phase 2 variables
	 */

	size_t *all_offset_source_file_to_send_phase2	= NULL;
	int *all_read_size_phase2_to_send				= NULL;
	size_t *all_offset_dest_file_phase2_to_send		= NULL;
	int *all_rank_phase2_to_send					= NULL;
	size_t *new_offset_dest_phase2 					= malloc(local_readNum* sizeof(size_t));
	size_t *new_offset_source_sorted_phase2 		= malloc(local_readNum* sizeof(size_t));
	int *new_read_size_phase2 						= malloc(local_readNum* sizeof(int));
	int *new_read_size_sorted_phase3 				= NULL;
	int *new_rank_phase2 							= malloc(local_readNum* sizeof(int));
	size_t *new_offset_dest_index_phase2 			= malloc(local_readNum* sizeof(size_t));
	size_t *pbs_local_source_offset					= NULL;
	size_t *pbs_local_source_offset_index			= NULL;

	size_t *all_offset_source_sorted_index_phase2	= NULL;
	//variables for the writing part

	//variables for MPI writes and read
	MPI_File out;
	char* path;
	double time_count;
	/*
	size_t *offset_source;
	offset_source = (size_t*)malloc(local_readNum*sizeof(size_t));
	offset_source[0] = 0;

	int *size_source;
	size_source = (int*)malloc(local_readNum*sizeof(int));
	size_source[0] = 0;
	*/

	/* TODO Change master_job at each iteration
	 * use master_rank_global_count
	 */

	char **data2;


	/* *****************************************************************************
	 * task: BEGIN PHASE 1
	 *
	 * In this phase we are going to sort the destination
	 * offset .and change the rank of the reads to tell them
	 * where to go to in order to be writen in the same block block
	 *
	 * Each job has new vector of offset read, of offset write
	 * and of read size :: new_offset_source, new_offset_dest,  new_read_size
	 *
	 * We need a new vector with the rank for sending the reads
	 * after reading.
	 *
	 * There is nothing really new compare with the sorting of the
	 * coordinates and same optimization could be done, except we compute
	 * new rank for each reads at the end
	 ******************************************************************************/


	if (rank == master_job_phase_2)
		fprintf(stderr, "Rank %d :::::[WRITE_ANY_DIM][PHASE 1] SORT OFFSET DESTINATION \n", rank);

	size_t *num_reads_per_jobs = (size_t *) malloc(total_num_proc * sizeof(size_t));

	MPI_Gather(
			&local_readNum,
			1,
			MPI_LONG_LONG_INT,
			&num_reads_per_jobs[rank - master_job_phase_2],
			1,
			MPI_LONG_LONG_INT,
			master_job_phase_1,
			COMM_WORLD
			);

		//vector of all offset in destination / source file / size of reads / ranks
		size_t *all_offset_dest_file_phase1		= NULL;
		size_t *all_offset_source_file_phase1	= NULL;
		int *all_read_size_phase1				= NULL;
		int *all_rank_phase1					= NULL;

		size_t total_num_read_phase1 = 0;
		MPI_Reduce(
				&local_readNum,
				&total_num_read_phase1,
				1,
				MPI_LONG_LONG_INT,
				MPI_SUM,
				master_job_phase_1,
				COMM_WORLD
				);

		if (rank == master_job_phase_1)
			assert(total_num_read_phase1 == total_num_read);

		if (rank == master_job_phase_1){
			all_offset_dest_file_phase1 	= malloc (total_num_read_phase1 * sizeof(size_t));
			all_offset_source_file_phase1 	= malloc (total_num_read_phase1 * sizeof(size_t));
			all_read_size_phase1 			= malloc (total_num_read_phase1 * sizeof(int));
			all_rank_phase1 				= malloc (total_num_read_phase1 * sizeof(int));

			size_t k =0;

			for (k = 0; k < total_num_read_phase1; k++){
				all_read_size_phase1[k] 			= 0;
				all_offset_dest_file_phase1[k] 		= 0;
				all_offset_source_file_phase1[k] 	= 0;
				all_rank_phase1[k] 					= 0;
			}
		}

		/*
		 * Phase 1: master_1 defines the following two vectors
		 */
		// vector of number of read per jobs
		size_t *num_reads_per_jobs_phase1 		= malloc(total_num_proc* sizeof(size_t));
		// vector of index that contains the cumulative sum of the number of reads
		size_t *start_num_reads_per_jobs_phase1 = malloc((total_num_proc + 1)*sizeof(size_t));

		/*
		 * Phase 2: master_2 receives all local_readNum and adds it to a local vector
		 */
		MPI_Gather(
				&local_readNum,
				1,
				MPI_LONG_LONG_INT,
				&num_reads_per_jobs_phase1[rank - master_job_phase_1],
				1,
				MPI_LONG_LONG_INT,
				master_job_phase_1,
				COMM_WORLD
				);

		if (rank == master_job_phase_1){

			start_num_reads_per_jobs_phase1[0] = 0;

			for (k = 1; k < (total_num_proc +1); k++){
				start_num_reads_per_jobs_phase1[k] = num_reads_per_jobs_phase1[k-1];
			}

			for (k = 1; k < total_num_proc; k++){
				size_t tmp 							= start_num_reads_per_jobs_phase1[k - 1];
				size_t tmp2 						= start_num_reads_per_jobs_phase1[k];
				start_num_reads_per_jobs_phase1[k] 	= tmp + tmp2;
			}
		}

		if (rank == master_job_phase_1){

			MPI_Status status;
			//we copy the first elements in
			int k=0;
			size_t st = start_num_reads_per_jobs_phase1[master_job_phase_1];

			for (k = 0; k < num_reads_per_jobs_phase1[master_job_phase_1]; k++){

				all_offset_dest_file_phase1[st] 	= new_offset_dest[k];
				all_offset_source_file_phase1[st] 	= new_offset_source[k];
				all_read_size_phase1[st] 			= new_read_size[k];
				all_rank_phase1[st] 				= new_rank[k];
				st++;
			}

			for(j = 0; j < total_num_proc; j++){

				if(j != master_job_phase_2){

					// first we care for ranks
					int *temp_buf 		= malloc(num_reads_per_jobs_phase1[j]* sizeof(int));
					int *temp_buf1 		= malloc(num_reads_per_jobs_phase1[j]* sizeof(int));
					size_t *temp_buf2 	= malloc(num_reads_per_jobs_phase1[j]* sizeof(size_t));
					size_t *temp_buf3 	= malloc(num_reads_per_jobs_phase1[j]* sizeof(size_t));

					MPI_Recv(temp_buf, num_reads_per_jobs_phase1[j], MPI_INT, j, 0, COMM_WORLD, &status);
					MPI_Recv(temp_buf1, num_reads_per_jobs_phase1[j], MPI_INT, j, 1, COMM_WORLD, &status);
					MPI_Recv(temp_buf2, num_reads_per_jobs_phase1[j], MPI_LONG_LONG_INT, j, 2, COMM_WORLD , &status);
					MPI_Recv(temp_buf3, num_reads_per_jobs_phase1[j], MPI_LONG_LONG_INT, j, 3, COMM_WORLD, &status);

					st=0;
					size_t st = start_num_reads_per_jobs_phase1[j];

					for (k = 0; k < num_reads_per_jobs_phase1[j]; k++){

						all_rank_phase1[st] 				= temp_buf[k];
						all_read_size_phase1[st] 			= temp_buf1[k];
						all_offset_source_file_phase1[st] 	= temp_buf2[k];
						all_offset_dest_file_phase1[st] 	= temp_buf3[k];
						st++;
					}

					free(temp_buf);
					free(temp_buf1);
					free(temp_buf2);
					free(temp_buf3);
				}
			}
		}
		else{
			MPI_Send(new_rank, local_readNum, MPI_INT, master_job_phase_1,  0,COMM_WORLD);
			MPI_Send(new_read_size, local_readNum, MPI_INT, master_job_phase_1,  1, COMM_WORLD);
			MPI_Send(new_offset_source, local_readNum, MPI_LONG_LONG_INT, master_job_phase_1,  2,COMM_WORLD);
			MPI_Send(new_offset_dest, local_readNum, MPI_LONG_LONG_INT, master_job_phase_1,  3, COMM_WORLD);
		}

		if (rank == master_job_phase_1){

			for ( j = 1; j < total_num_read_phase1; j++){
				assert ( all_offset_dest_file_phase1[j] != 0 );
				assert ( all_offset_source_file_phase1[j] != 0 );
				assert ( all_read_size_phase1[j] != 0 );
				assert ( all_rank_phase1[j] <= total_num_proc);
			}
		}

		/**************************/
		// We free some variable
		/**************************/

		free(new_read_size);
		free(new_offset_dest);
		free(new_offset_source);
		free(new_rank);

		/*
		 * In this section we implement a parallel Bitonic sort
		 * algorithm.
		 * Input are
		 * all_read_size_phase1,
		 * all_offset_dest_index_phase1
		 * all_offset_dest_file_phase1
		 *
		 */


		// the master rank compute the number of
		// dimension is the number of processors where we
		// perform the bitonic sort
		// int dimensions = (int)(log2(num_processes));
		// find next ( must be greater) power, and go one back

		// we broadcast the total number of reads to each rank
		MPI_Bcast(&total_num_read_phase1, 1, MPI_LONG_LONG_INT, master_job_phase_1, COMM_WORLD );

		assert(total_num_read_phase1 == total_num_read);
		/*
		int dimensions = 1;
		while (dimensions <= total_num_proc)
			dimensions <<= 1;

		dimensions >>= 1;
		*/
		// the master rank compute the number of
		// dimension is the number of processors where we
		// perform the bitonic sort
		// int dimensions = (int)(log2(num_processes));
		// find next ( must be greater) power, and go one back

		// we broadcast the total number of reads to each rank
		// we compute the length of the vector to recieve
		// the we split all_offset_dest_index_phase1 among the dimension processors
		size_t pbs_num_offsets_to_recieve_phase1 = total_num_read_phase1/dimensions;

		// for the last processor will recieve the rest of the division
		size_t pbs_num_offsets_to_recieve_left_phase1  =  total_num_read_phase1 - pbs_num_offsets_to_recieve_phase1*dimensions;

		//fprintf(stderr, "pbs_num_offsets_to_recieve_left =  %zu \n", pbs_num_offsets_to_recieve_left);

		/*
		 * Here we reduce the dimension.
		 * We do that because we don't want a job to finish with zero vector
		 * after the bitonic.
		 *
		 * In the future we shall use a bruck to dispatch the read between dimensions
		 */
		while((pbs_num_offsets_to_recieve_left_phase1 * dimensions) > pbs_num_offsets_to_recieve_phase1){
			dimensions >>= 1;
			pbs_num_offsets_to_recieve_phase1 		= total_num_read_phase1/dimensions;
			pbs_num_offsets_to_recieve_left_phase1  = total_num_read_phase1 - pbs_num_offsets_to_recieve_phase1*dimensions;
		}

		// we compute a vector of size dimensions which contain the number
		// of reads to send
		size_t *pbs_local_num_read_per_job_phase1 = malloc(dimensions * sizeof(size_t));

		// each job create a vector with with the length
		// of the offset vector for each job in [0, dimension[
		for (j = 0; j < dimensions ; j++){
			// we add pbs_num_offsets_to_recieve_left
			// because we want the vector to have the same size
			pbs_local_num_read_per_job_phase1[j] = pbs_num_offsets_to_recieve_phase1 + pbs_num_offsets_to_recieve_left_phase1;
		}
		size_t *pbs_start_num_offset_per_jobs_phase1 = malloc((dimensions + 1) * sizeof(size_t));

		// the master job compute the start index of the element
		// to dispatch
		if (rank == master_job_phase_1){

			pbs_start_num_offset_per_jobs_phase1[0] = 0;

			for (k = 1; k < (dimensions +1); k++){
				pbs_start_num_offset_per_jobs_phase1[k] = pbs_local_num_read_per_job_phase1[k-1];
			}
			for (k = 1; k < dimensions; k++){
				size_t tmp 	= pbs_start_num_offset_per_jobs_phase1[k - 1];
				size_t tmp2 = pbs_start_num_offset_per_jobs_phase1[k];
				// we remove the left over reads
				pbs_start_num_offset_per_jobs_phase1[k] = tmp + tmp2 - pbs_num_offsets_to_recieve_left_phase1;
			}
		}

		// the processors master_job_phase_1 send the output offset
		// to all the rank in [0-dimension]
		if (rank < dimensions){

			// pbs_local_offset_to_sort is a table containing the unsorted
			// destination offset
			pbs_local_dest_offset = malloc(sizeof(size_t) * pbs_local_num_read_per_job_phase1[rank]);
			//now the master send

			if ( rank != master_job_phase_1 ){
					MPI_Recv(
							pbs_local_dest_offset,
							pbs_local_num_read_per_job_phase1[rank],
							MPI_LONG_LONG_INT,
							master_job_phase_1,
							0,
							COMM_WORLD, &status
							);

			}
			else {
				//first we copy the data from the master job
				size_t ind = pbs_start_num_offset_per_jobs_phase1[master_job_phase_1];

				for (k = 0; k < pbs_local_num_read_per_job_phase1[master_job_phase_1]; k++){
					pbs_local_dest_offset[k] = all_offset_dest_file_phase1[ind];
					ind++;
				}

				for(j = 0; j < dimensions; j++){
					if (j != master_job_phase_1){
						MPI_Send(
								&all_offset_dest_file_phase1[pbs_start_num_offset_per_jobs_phase1[j]],
								pbs_local_num_read_per_job_phase1[j],
								MPI_LONG_LONG_INT,
								j,
								0,
								COMM_WORLD
								);
					}
				}
			}

			// we build pbs_local_dest_offset_index
			pbs_local_dest_offset_index = malloc(pbs_local_num_read_per_job_phase1[rank]*sizeof(size_t));

			for (j = 0; j < pbs_local_num_read_per_job_phase1[rank]; j++){

				if (rank == master_job_phase_1){
					pbs_local_dest_offset_index[j] = j + pbs_local_num_read_per_job_phase1[rank]*rank;
				}
				else{
					pbs_local_dest_offset_index[j] = j + pbs_local_num_read_per_job_phase1[rank]*rank -
							(rank*pbs_num_offsets_to_recieve_left_phase1);
				}
			}

			for ( j = 0; j < pbs_local_num_read_per_job_phase1[rank]; j++){
				assert ( pbs_local_dest_offset[j] != 0 );
			}


			// now each rank from [0, dimension[
			// is going to bitonic sort
			// input are:
			// pbs_local_dest_offset
			// pbs_local_dest_offset_index

			// we call the parallel bitonic sort

			ParallelBitonicSort(
					COMM_WORLD,
					rank,
					dimensions,
					pbs_local_dest_offset,
					pbs_local_dest_offset_index,
					pbs_local_num_read_per_job_phase1[rank],
					pbs_num_offsets_to_recieve_left_phase1);


			for(j = 0; j < pbs_local_num_read_per_job_phase1[rank]; j++){
				assert(pbs_local_dest_offset_index[j] <= total_num_read_phase1);
			}

			time_count = MPI_Wtime();

			//we compute a new total number of reads
			size_t total_num_read_after_bitonic_sort = 0;
			for (k = 0; k < dimensions; k++){
				total_num_read_after_bitonic_sort += pbs_local_num_read_per_job_phase1[k];
			}

			// now we gather all the pbs_local_dest_offset_index
			// and pbs_local_dest_offset in 2 vectors
			// all_offset_dest_sorted_phase1
			// all_offset_index_phase_1

			//we allocate vector to send
			// we remove zero
			size_t start_index=0;
			while (pbs_local_dest_offset[start_index] == 0)
				start_index++;


			pbs_local_num_read_per_job_phase1[rank] -=  start_index;

			all_offset_dest_file_to_send_phase1 = (size_t *)malloc(sizeof(size_t) * total_num_read_phase1);
			all_offset_dest_sorted_index_phase1 = (size_t *)malloc(sizeof(size_t) * total_num_read_phase1);

			if (rank == master_job_phase_1){

				pbs_start_num_offset_per_jobs_phase1[0] = 0;

				for (k = 1; k < (dimensions +1); k++){
					pbs_start_num_offset_per_jobs_phase1[k] = pbs_local_num_read_per_job_phase1[k-1];
				}
				for (k = 1; k < dimensions; k++){
					size_t tmp = pbs_start_num_offset_per_jobs_phase1[k - 1];
					size_t tmp2 = pbs_start_num_offset_per_jobs_phase1[k];
					// we remove the left over reads
					pbs_start_num_offset_per_jobs_phase1[k] = tmp + tmp2;
				}
			}

			time_count = MPI_Wtime();

			// we gather the offset dest sorted
			chosen_split_rank_gather_size_t(
					COMM_WORLD,
					rank,
					dimensions,
					master_job_phase_1,
					pbs_local_num_read_per_job_phase1[rank],
					pbs_local_num_read_per_job_phase1,
					pbs_start_num_offset_per_jobs_phase1,
					all_offset_dest_file_to_send_phase1,
					pbs_local_dest_offset,
					start_index);

			if (rank == master_job_phase_1){
				for ( j = 1; j < total_num_read_phase1; j++){
					assert ( all_offset_dest_file_to_send_phase1[j] != 0 );
				}
				for ( j = 0; j < total_num_read_phase1-1; j++){
					assert( all_offset_dest_file_to_send_phase1[j] < all_offset_dest_file_to_send_phase1[j + 1] );
				}
			}

			chosen_split_rank_gather_size_t(
							COMM_WORLD,
							rank,
							dimensions,
							master_job_phase_1,
							pbs_local_num_read_per_job_phase1[rank],
							pbs_local_num_read_per_job_phase1,
							pbs_start_num_offset_per_jobs_phase1,
							all_offset_dest_sorted_index_phase1,
							pbs_local_dest_offset_index,
							start_index);


			if (rank == master_job_phase_1){

				// now we apply the new index to all the
				all_read_size_phase1_to_send 			= malloc(total_num_read_phase1*sizeof(int));
				all_offset_source_file_phase1_to_send   = malloc(total_num_read_phase1*sizeof(size_t));
				all_rank_phase1_to_send 				= malloc(total_num_read_phase1*sizeof(int));

				for(j = 0; j < total_num_read_phase1; j++){

					all_offset_source_file_phase1_to_send[j] = all_offset_source_file_phase1[all_offset_dest_sorted_index_phase1[j]];
					all_read_size_phase1_to_send[j] 		 = all_read_size_phase1[all_offset_dest_sorted_index_phase1[j]];

				}

				//now we change the rank
				// we initialize all_offset_rank_to_send constains
				// the rank of the sorted read
				size_t total = 0;
				for(j = 0; j < total_num_proc; j++){

					total += num_reads_per_jobs_phase1[j];
					for (k = 0; k < num_reads_per_jobs_phase1[j]; k++){
						all_rank_phase1_to_send[start_num_reads_per_jobs_phase1[j] + k] = j;
						}
				}
			} // end if (rank == master_job_phase_1)

		} //end if (rank < dimensions)

		MPI_Barrier(COMM_WORLD);

		if (rank < dimensions){

			free(pbs_local_dest_offset);
			free(pbs_local_dest_offset_index);
			free(all_offset_dest_sorted_index_phase1);


		}
		free(pbs_local_num_read_per_job_phase1);
		free(pbs_start_num_offset_per_jobs_phase1);

		if (rank == master_job_phase_1){

			free(all_offset_dest_file_phase1);
			free(all_read_size_phase1);
			free(all_rank_phase1);
			free(all_offset_source_file_phase1);
		}

		 //  task Phase 1: Dispatch everything

		if (rank != master_job_phase_1){

			MPI_Recv( new_offset_dest_phase1, local_readNum, MPI_LONG_LONG_INT, master_job_phase_1, 0, COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv( new_offset_source_phase1, local_readNum, MPI_LONG_LONG_INT, master_job_phase_1, 1,COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv( new_read_size_phase1, local_readNum, MPI_INT, master_job_phase_1, 2, COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv( new_rank_phase1, local_readNum, MPI_INT, master_job_phase_1, 3, COMM_WORLD, MPI_STATUS_IGNORE);

		}
		else {
			size_t k=0;
			size_t ind = start_num_reads_per_jobs_phase1[master_job_phase_1];


			for (k = 0; k < (num_reads_per_jobs_phase1[master_job_phase_1]); k++){

				new_offset_dest_phase1[k] 	= all_offset_dest_file_to_send_phase1[ind];
				new_offset_source_phase1[k] = all_offset_source_file_phase1_to_send[ind];
				new_read_size_phase1[k] 	= all_read_size_phase1_to_send[ind];
				new_rank_phase1[k] 			= all_rank_phase1_to_send[ind];
				ind++;
			}

			for(j = 0; j < total_num_proc; j++){

				if (j != master_job_phase_1){

					MPI_Send(
							&all_offset_dest_file_to_send_phase1[start_num_reads_per_jobs_phase1[j]],
							num_reads_per_jobs_phase1[j],
							MPI_LONG_LONG_INT,
							j,
							0,
							COMM_WORLD
							);

					MPI_Send(
							&all_offset_source_file_phase1_to_send[start_num_reads_per_jobs_phase1[j]],
							num_reads_per_jobs_phase1[j],
							MPI_LONG_LONG_INT,
							j,
							1,
							COMM_WORLD
							);

					MPI_Send(
							&all_read_size_phase1_to_send[start_num_reads_per_jobs_phase1[j]],
							num_reads_per_jobs_phase1[j],
							MPI_INT,
							j,
							2,
							COMM_WORLD
							);

					MPI_Send(
							&all_rank_phase1_to_send[start_num_reads_per_jobs_phase1[j]],
							num_reads_per_jobs_phase1[j],
							MPI_INT,
							j,
							3,
							COMM_WORLD
							);

				}
			}
		}
		/*
		 * Phase 1: Clean memory
		 */

		if (rank < dimensions){
			free(all_offset_dest_file_to_send_phase1);
		}

		if (rank == master_job_phase_1){
			free(all_read_size_phase1_to_send);
			free(all_rank_phase1_to_send);
			free(all_offset_source_file_phase1_to_send);
		}

		free(start_num_reads_per_jobs_phase1);
		free(num_reads_per_jobs_phase1);


	/* *****************************************************************************
	 * task: BEGIN PHASE 2
	 *
	 * In this phase we are going to sort the source
	 * offset . In order to read consecutive blocks each
	 * jobs.
	 *
	 * Each job has new vector of offset read, of offset write
	 * and of read size :: new_offset_source, new_offset_dest,  new_read_size
	 *
	 * We need a new vector with the rank for sending the reads
	 * after reading.
	 *
	 ******************************************************************************/

	//size_t *num_reads_per_jobs = (size_t *) malloc(num_proc * sizeof(size_t));
	MPI_Gather(
			&local_readNum,
			1,
			MPI_LONG_LONG_INT,
			&num_reads_per_jobs[rank - master_job_phase_2],
			1,
			MPI_LONG_LONG_INT,
			master_job_phase_2,
			COMM_WORLD
			);

	if (rank == master_job_phase_2)
		fprintf(stderr, "Rank %d :::::[WRITE_ANY_DIM][PHASE 2] SORT OFFSET SOURCES \n", rank);


	//vector of all offset in destination / source file / size of reads / ranks
	size_t *all_offset_dest_file_phase2		= NULL;
	size_t *all_offset_source_file_phase2	= NULL;
	int *all_read_size_phase2				= NULL;
	int *all_rank_phase2					= NULL;

	/*
	 * Phase 2: master_2 gets the total of reads
	 * Improvement: see if the total num read is the same of the phase1
	 */


	size_t total_num_read_phase2 = 0;
	MPI_Reduce(&local_readNum, &total_num_read_phase2, 1, MPI_LONG_LONG_INT, MPI_SUM, master_job_phase_2, COMM_WORLD);

	if (rank == master_job_phase_2)
		assert(total_num_read_phase2 == total_num_read);

	if (rank == master_job_phase_2){
		all_offset_dest_file_phase2 	= malloc (total_num_read_phase2 * sizeof(size_t));
		all_offset_source_file_phase2 	= malloc (total_num_read_phase2 * sizeof(size_t));
		all_read_size_phase2 			= malloc (total_num_read_phase2 * sizeof(int));
		all_rank_phase2 				= malloc (total_num_read_phase2 * sizeof(int));

		size_t k =0;

		for (k = 0; k < total_num_read_phase2; k++){
			all_read_size_phase2[k] 			= 0;
			all_offset_dest_file_phase2[k] 		= 0;
			all_offset_source_file_phase2[k] 	= 0;
			all_rank_phase2[k] 					= 0;
		}
	}

	// vector of number of read per jobs
	size_t *num_reads_per_jobs_phase2 = malloc(total_num_proc* sizeof(size_t));
	// vector of index that contains the cumulative sum of the number of reads
	size_t *start_num_reads_per_jobs_phase2 = malloc((total_num_proc + 1)*sizeof(size_t));

	/*
	 * Phase 2: master_2 receives all local_readNum and adds it to a local vector
	 */
	MPI_Gather(
			&local_readNum,
			1,
			MPI_LONG_LONG_INT,
			&num_reads_per_jobs_phase2[rank - master_job_phase_2],
			1,
			MPI_LONG_LONG_INT,
			master_job_phase_2,
			COMM_WORLD
			);


	if (rank == master_job_phase_2){

		start_num_reads_per_jobs_phase2[0] = 0;

		for (k = 1; k < (total_num_proc +1); k++){
			start_num_reads_per_jobs_phase2[k] = num_reads_per_jobs_phase2[k-1];
		}

		for (k = 1; k < total_num_proc; k++){
			size_t tmp 	= start_num_reads_per_jobs_phase2[k - 1];
			size_t tmp2 = start_num_reads_per_jobs_phase2[k];
			start_num_reads_per_jobs_phase2[k] = tmp + tmp2;
		}
	}

	/*
	 * split_chosen_rank
	 * Collect sizes, coordinates and offsets
	 * in all_vector
	 */

	if (rank == master_job_phase_2){

		MPI_Status status;
		//we copy the first elements in
		int k=0;
		size_t st = start_num_reads_per_jobs_phase2[master_job_phase_2];

		for (k = 0; k < num_reads_per_jobs[master_job_phase_2]; k++){

			all_offset_dest_file_phase2[st] 	= new_offset_dest_phase1[k];
			all_offset_source_file_phase2[st] 	= new_offset_source_phase1[k];
			all_read_size_phase2[st] 			= new_read_size_phase1[k];
			all_rank_phase2[st] 				= new_rank_phase1[k];
			st++;
		}

		for(j = 0; j < total_num_proc; j++){

			if(j != master_job_phase_2){

				// first we care for ranks
				int *temp_buf 		= malloc(num_reads_per_jobs_phase2[j]* sizeof(int));
				int *temp_buf1 		= malloc(num_reads_per_jobs_phase2[j]* sizeof(int));
				size_t *temp_buf2 	= malloc(num_reads_per_jobs_phase2[j]* sizeof(size_t));
				size_t *temp_buf3 	= malloc(num_reads_per_jobs_phase2[j]* sizeof(size_t));

				MPI_Recv(temp_buf, num_reads_per_jobs_phase2[j], MPI_INT, j, 0, COMM_WORLD, &status);
				MPI_Recv(temp_buf1, num_reads_per_jobs_phase2[j], MPI_INT, j, 1, COMM_WORLD, &status);
				MPI_Recv(temp_buf2, num_reads_per_jobs_phase2[j], MPI_LONG_LONG_INT, j, 2, COMM_WORLD , &status);
				MPI_Recv(temp_buf3, num_reads_per_jobs_phase2[j], MPI_LONG_LONG_INT, j, 3, COMM_WORLD, &status);

				st=0;
				size_t st = start_num_reads_per_jobs_phase2[j];

				for (k = 0; k < num_reads_per_jobs_phase2[j]; k++){

					all_rank_phase2[st] 				= temp_buf[k];
					all_read_size_phase2[st] 			= temp_buf1[k];
					all_offset_source_file_phase2[st] 	= temp_buf2[k];
					all_offset_dest_file_phase2[st] 	= temp_buf3[k];
					st++;
				}

				free(temp_buf);
				free(temp_buf1);
				free(temp_buf2);
				free(temp_buf3);
			}
		}
	}
	else{
		MPI_Send(new_rank_phase1, local_readNum, MPI_INT, master_job_phase_2,  0,COMM_WORLD);
		MPI_Send(new_read_size_phase1, local_readNum, MPI_INT, master_job_phase_2,  1, COMM_WORLD);
		MPI_Send(new_offset_source_phase1, local_readNum, MPI_LONG_LONG_INT, master_job_phase_2,  2,COMM_WORLD);
		MPI_Send(new_offset_dest_phase1, local_readNum, MPI_LONG_LONG_INT, master_job_phase_2,  3, COMM_WORLD);
	}

	if (rank == master_job_phase_2){
		for ( j = 1; j < total_num_read_phase2; j++){
			assert ( all_offset_dest_file_phase2[j] 	!= 0 );
			assert ( all_offset_source_file_phase2[j] 	!= 0 );
			assert ( all_read_size_phase2[j] 			!= 0 );
			assert ( all_rank_phase2[j] 				<= total_num_proc);
		}
	}

	/**************************/
	// We free some variable
	/**************************/

	free(new_read_size_phase1);
	free(new_offset_dest_phase1);
	free(new_offset_source_phase1);
	free(new_rank_phase1);
	free(num_reads_per_jobs_phase2);

	/*
	 * task phase 3: bitonic sort of source offset
	 * In this section we implement a parallel Bitonic sort
	 * algorithm.
	 * Input are
	 * all_read_size_phase1,
	 * all_offset_dest_index_phase1
	 * all_offset_dest_file_phase1
	 *
	 */

	// the master rank compute the number of
	// dimension is the number of processors where we
	// perform the bitonic sort
	// int dimensions = (int)(log2(num_processes));
	// find next ( must be greater) power, and go one back

	// we broadcast the total number of reads to each rank
	MPI_Bcast(&total_num_read_phase2, 1, MPI_LONG_LONG_INT, master_job_phase_2, COMM_WORLD );

	assert(total_num_read_phase2 == total_num_read);

	// we broadcast the total number of reads to each rank
	// we compute the length of the vector to recieve
	// the we split all_offset_dest_index_phase1 among the dimension processors
	size_t pbs_num_offsets_to_recieve = total_num_read_phase2/dimensions;
	// for the last processor will recieve the rest of the division
	size_t pbs_num_offsets_to_recieve_left  =  total_num_read_phase2 - pbs_num_offsets_to_recieve*dimensions;
	//fprintf(stderr, "pbs_num_offsets_to_recieve_left =  %zu \n", pbs_num_offsets_to_recieve_left);

	while((pbs_num_offsets_to_recieve_left * dimensions) > pbs_num_offsets_to_recieve){
		dimensions >>= 1;
		pbs_num_offsets_to_recieve = total_num_read_phase2/dimensions;
		pbs_num_offsets_to_recieve_left  = total_num_read_phase2 - pbs_num_offsets_to_recieve*dimensions;
	}

	// we compute a vector of size dimensions which contain the number
	// of reads to send
	size_t *pbs_local_num_read_per_job_phase2 = malloc(dimensions * sizeof(size_t));

	// each job create a vector with with the length
	// of the offset vector for each job in [0, dimension[
	for (j = 0; j < dimensions ; j++){
		// we add pbs_num_offsets_to_recieve_left
		// because we want the vector to have the same size
		pbs_local_num_read_per_job_phase2[j] = pbs_num_offsets_to_recieve + pbs_num_offsets_to_recieve_left;
	}


	size_t *pbs_start_num_offset_per_jobs_phase2 = malloc((dimensions + 1) * sizeof(size_t));

	// the master job compute the start index of the element
	// to dispatch
	if (rank == master_job_phase_2){

		pbs_start_num_offset_per_jobs_phase2[0] = 0;

		for (k = 1; k < (dimensions +1); k++){
			pbs_start_num_offset_per_jobs_phase2[k] = pbs_local_num_read_per_job_phase2[k-1];
		}
		for (k = 1; k < dimensions; k++){
			size_t tmp = pbs_start_num_offset_per_jobs_phase2[k - 1];
			size_t tmp2 = pbs_start_num_offset_per_jobs_phase2[k];
			// we remove the left over reads
			pbs_start_num_offset_per_jobs_phase2[k] = tmp + tmp2 - pbs_num_offsets_to_recieve_left;
		}
	}

	// the processors master_job_phase_2 send the source offset
	// to all the rank in [0-dimension]
	if (rank < dimensions){

		// pbs_local_offset_to_sort is a table containing the unsorted
		// destination offset
		pbs_local_source_offset = malloc(sizeof(size_t) * pbs_local_num_read_per_job_phase2[rank]);
		//now the master send

		if ( rank != master_job_phase_2 ){
				MPI_Recv(
						pbs_local_source_offset,
						pbs_local_num_read_per_job_phase2[rank], MPI_LONG_LONG_INT,
						master_job_phase_2,
						0,
						COMM_WORLD,
						&status
						);

		}
		else {
			//first we copy the data from the master job
			size_t ind = pbs_start_num_offset_per_jobs_phase2[master_job_phase_2];

			for (k = 0; k < pbs_local_num_read_per_job_phase2[master_job_phase_2]; k++){
				pbs_local_source_offset[k] = all_offset_source_file_phase2[ind];
				ind++;
			}

			for(j = 0; j < dimensions; j++){
				if (j != master_job_phase_2){
					MPI_Send(
							&all_offset_source_file_phase2[pbs_start_num_offset_per_jobs_phase2[j]],
							pbs_local_num_read_per_job_phase2[j],
							MPI_LONG_LONG_INT,
							j,
							0,
							COMM_WORLD
							);
				}
			}
		}
		// we build pbs_local_dest_offset_index
		pbs_local_source_offset_index = malloc(pbs_local_num_read_per_job_phase2[rank]*sizeof(size_t));

		//fprintf(stderr, "Rank %d :::::[WRITE] pbs_num_offsets_to_recieve_left = %zu \n", rank, pbs_num_offsets_to_recieve_left);

		for (j = 0; j < pbs_local_num_read_per_job_phase2[rank]; j++){

			if (rank == master_job_phase_2){
				pbs_local_source_offset_index[j] = j + pbs_local_num_read_per_job_phase2[rank]*rank;
			}
			else{
				pbs_local_source_offset_index[j] = j + pbs_local_num_read_per_job_phase2[rank]*rank - (rank*pbs_num_offsets_to_recieve_left);
			}
		}

		for ( j = 0; j < pbs_local_num_read_per_job_phase2[rank]; j++){
			assert ( pbs_local_source_offset[j] != 0 );
		}


		// now each rank from [0, dimension[
		// is going to bitonic sort
		// input are:
		// pbs_local_dest_offset
		// pbs_local_dest_offset_index

		// we call the parallel bitonic sort
		time_count = MPI_Wtime();
		ParallelBitonicSort(
				COMM_WORLD,
				rank,
				dimensions,
				pbs_local_source_offset,
				pbs_local_source_offset_index,
				pbs_local_num_read_per_job_phase2[rank],
				pbs_num_offsets_to_recieve_left);

		for(j = 1; j < pbs_local_num_read_per_job_phase2[rank]; j++){
			assert(pbs_local_source_offset[j-1] <= pbs_local_source_offset[j]);
			assert(pbs_local_source_offset_index[j] <= total_num_read_phase2);
		}

		//we compute a new total number of reads
		size_t total_num_read_after_bitonic_sort = 0;
		for (k = 0; k < dimensions; k++){
			total_num_read_after_bitonic_sort += pbs_local_num_read_per_job_phase2[k];
		}

		// now we gather all the pbs_local_dest_offset_index
		// and pbs_local_dest_offset in 2 vectors
		// all_offset_dest_sorted_phase1
		// all_offset_index_phase_1

		//we allocate vector to send
		// we remove zero
		size_t start_index=0;
		while (pbs_local_source_offset[start_index] == 0)
			start_index++;

		pbs_local_num_read_per_job_phase2[rank] -=  start_index;

		all_offset_source_file_to_send_phase2 = malloc(sizeof(size_t) * total_num_read_phase2);
		all_offset_source_sorted_index_phase2 = malloc(sizeof(size_t) * total_num_read_phase2);

		if (rank == master_job_phase_2){

			pbs_start_num_offset_per_jobs_phase2[0] = 0;

			for (k = 1; k < (dimensions +1); k++){
				pbs_start_num_offset_per_jobs_phase2[k] = pbs_local_num_read_per_job_phase2[k-1];
			}
			for (k = 1; k < dimensions; k++){
				size_t tmp = pbs_start_num_offset_per_jobs_phase2[k - 1];
				size_t tmp2 = pbs_start_num_offset_per_jobs_phase2[k];
				// we remove the left over reads
				pbs_start_num_offset_per_jobs_phase2[k] = tmp + tmp2;
			}
		}

		time_count = MPI_Wtime();
		// we gather the offset source sorted
		chosen_split_rank_gather_size_t(
				COMM_WORLD,
				rank,
				dimensions,
				master_job_phase_2,
				pbs_local_num_read_per_job_phase2[rank],
				pbs_local_num_read_per_job_phase2,
				pbs_start_num_offset_per_jobs_phase2,
				all_offset_source_file_to_send_phase2,
				pbs_local_source_offset,
				start_index);

		if (rank == master_job_phase_2){
			for ( j = 1; j < total_num_read_phase2; j++){
				assert ( all_offset_source_file_to_send_phase2[j] != 0 );
			}
			for ( j = 0; j < total_num_read_phase2-1; j++){
				assert( all_offset_source_file_to_send_phase2[j] < all_offset_source_file_to_send_phase2[j + 1] );
			}
		}

		chosen_split_rank_gather_size_t(
						COMM_WORLD,
						rank,
						dimensions,
						master_job_phase_2,
						pbs_local_num_read_per_job_phase2[rank],
						pbs_local_num_read_per_job_phase2,
						pbs_start_num_offset_per_jobs_phase2,
						all_offset_source_sorted_index_phase2,
						pbs_local_source_offset_index,
						start_index
						);

		if (rank == master_job_phase_2){
			fprintf(stderr, "Rank %d :::::[WRITE_ANY_DIM][PHASE 3] Time for gathering source offset and index %f\n", rank, MPI_Wtime() - time_count);


			// now we apply the new index to all the
			all_read_size_phase2_to_send 		= malloc(total_num_read_phase2*sizeof(int));
			all_offset_dest_file_phase2_to_send = malloc(total_num_read_phase2*sizeof(size_t));
			all_rank_phase2_to_send 			= malloc(total_num_read_phase2*sizeof(int));

			for(j = 0; j < total_num_read_phase2; j++){

				all_offset_dest_file_phase2_to_send[j] = all_offset_dest_file_phase2[all_offset_source_sorted_index_phase2[j]];
				all_read_size_phase2_to_send[j] 	   = all_read_size_phase2[all_offset_source_sorted_index_phase2[j]];
				all_rank_phase2_to_send[j] 			   = all_rank_phase2[all_offset_source_sorted_index_phase2[j]];
			}

			if (rank == master_job_phase_2){
				for(j = 0; j < total_num_read_phase2; j++){
					assert ( all_offset_dest_file_phase2_to_send[j]	!= 0 );
					assert ( all_read_size_phase2_to_send[j] 		!= 0 );
					assert ( all_rank_phase2_to_send[j] 			<= total_num_proc );
				}
			}

			size_t total = 0;
			for(j = 0; j < total_num_proc; j++){
				total += num_reads_per_jobs[j];
			}
			assert(total == total_num_read_phase2);

		} // end if (rank == master_job_phase_2)

	} //end if (rank < dimensions)

	MPI_Barrier(COMM_WORLD);

	if (rank < dimensions){
		free(pbs_local_source_offset);
		free(pbs_local_source_offset_index);
		free(all_offset_source_sorted_index_phase2);
	}

	free(pbs_local_num_read_per_job_phase2);
	free(pbs_start_num_offset_per_jobs_phase2);

	if (rank == master_job_phase_2){
		free(all_offset_dest_file_phase2);
		free(all_read_size_phase2);
		free(all_rank_phase2);
		free(all_offset_source_file_phase2);
	}

	//  task Phase 2: Dispatch everything


	//Phase 2: Send all_offsets_dest from master_2 to all
	send_size_t_master_to_all(rank, total_num_proc, master_job_phase_2, local_readNum, num_reads_per_jobs, start_num_reads_per_jobs_phase2,
			all_offset_dest_file_phase2_to_send, new_offset_dest_phase2);

	//Phase 2: Send all_offsets_source from master_2 to all
	send_size_t_master_to_all(rank, total_num_proc, master_job_phase_2, local_readNum, num_reads_per_jobs, start_num_reads_per_jobs_phase2,
			all_offset_source_file_to_send_phase2, new_offset_source_sorted_phase2);

	//Phase 2: Send all_size from master_2 to all
	send_int_master_to_all(rank, total_num_proc, master_job_phase_2, local_readNum, num_reads_per_jobs, start_num_reads_per_jobs_phase2,
			all_read_size_phase2_to_send, new_read_size_phase2);

	//Phase 2: Send all_rank from master_2 to all
	send_int_master_to_all(rank, total_num_proc, master_job_phase_2, local_readNum, num_reads_per_jobs, start_num_reads_per_jobs_phase2,
			all_rank_phase2_to_send, new_rank_phase2);

	/*
	 * FOR DEBUG
	 *
	 *
	for ( j = 0; j < (local_readNum - 1); j++){
		assert( new_offset_source_sorted_phase2[j] < new_offset_source_sorted_phase2[j + 1]);
	}

	for ( j = 0; j < local_readNum ; j++){
		assert(new_offset_dest_phase2[j] != 0);
		assert(new_read_size_phase2[j] != 0);
		assert(new_rank_phase2[j] < num_proc);

	}
	 */

	/*
	 * Phase 2: Clean memory
	 */
	if (rank == master_job_phase_2){
		//we free pointers
		free(all_offset_dest_file_phase2_to_send);
		free(all_read_size_phase2_to_send);
		free(all_rank_phase2_to_send);
	}

	if (rank < dimensions)
		free(all_offset_source_file_to_send_phase2);

	free(start_num_reads_per_jobs_phase2);
	free(num_reads_per_jobs);

	/***************************************************/
	/*
	 * In this part we are going to read
	 * the input reads according to sorting offset
	 *
	 *	new_offset_dest_phase2 (offset in the destination file) not sorted
	 *	new_rank_phase2 (rank per read) not sorted
	 *	new_read_size_phase2 (size of read ) not sorted
	 *	new_offset_source_sorted_phase2 (offset in the source file of the reads) sorted
	 *
	 */
	/***************************************************/


	/* *****************************************
	 * task Reading phase
	 * we compute the size of
	 * data for each rank and we put it in
	 * buffs and buffs_by_proc is
	 * the size of the buffer to send
	 *******************************************/
	size_t *number_of_reads_by_procs = (size_t*)calloc( total_num_proc, sizeof(size_t));
	size_t *buffs_by_procs 			 = (size_t*)calloc( total_num_proc, sizeof(size_t));

	time_count = MPI_Wtime();
	int num_proc = total_num_proc;

	/*
	 * first we pack data in data_pack
	 * to have consecutive reads of the
	 * chromosoms
	 */


	time_count = MPI_Wtime();

	int m;
	size_t i;
	size_t new_data_sz = 0;


	char *data_pack;
	char *tmp_tab;
	//we compute the size of data_pack
	for (k = 0; k < local_readNum; k++){
		new_data_sz += new_read_size_phase2[k];
	}

	data_pack = malloc(new_data_sz +1);
	data_pack[new_data_sz] = 0;

	char *q = data;
	char *p = data_pack;
	size_t offset_in_data = 0;
	int pos = 0;
	//we compute the new offset of reads in data buffer
	//we remove the start offset in the file

	//we copy elements from data in data_pack
	for (k=0; k < local_readNum; k++){
		pos = 0;
		offset_in_data = new_offset_source_sorted_phase2[k] - start_offset_in_file;
		q = data + offset_in_data;
		while (*q && (pos < new_read_size_phase2[k])) {	*p=*q; q++;	p++; pos++;}
	}

	int res;

	/*
	 * We unpack in a loop the same way
	 */

	MPI_Datatype dt_data;

	//The data in which what we read will be kept
	data2 = (char**)malloc( num_proc * sizeof(char*));

	// we compute the size of
	// data for each rank and we put it in
	// buffs and buffs_by_proc is
	// the size of the buffer to send
	size_t *buffs_by_procs2 = (size_t*)calloc( num_proc, sizeof(size_t));

	for(m = 0; m < local_readNum; m++)
	{
		buffs_by_procs2[new_rank_phase2[m]] += new_read_size_phase2[m];
		number_of_reads_by_procs[new_rank_phase2[m]]++;
	}

	for(m = 0; m < num_proc; m++)
	{
		buffs_by_procs[(rank + m)%num_proc] = buffs_by_procs2[(rank - m + num_proc)%num_proc];
	}

	free(buffs_by_procs2);

	//Allocate data and initialization
	for(m = 0; m < num_proc; m++)
	{
		data2[m] = (char*)malloc(buffs_by_procs[m]*sizeof(char) + 1);
		data2[m][buffs_by_procs[m]] = 0;
	}

	//Variable for datatype struct
	MPI_Aint *indices = (MPI_Aint *)malloc(local_readNum*sizeof(MPI_Aint));
	int *blocklens = (int *)malloc(local_readNum*sizeof(int));

	MPI_Datatype *oldtypes = (MPI_Datatype *)malloc(local_readNum*sizeof(MPI_Datatype));
	MPI_Aint adress_to_write_in_data_by_element[num_proc];
	for(i = 0; i < num_proc; i++){
		MPI_Get_address(data2[(rank-i+num_proc)%num_proc], &adress_to_write_in_data_by_element[(rank+i)%num_proc]);
	}
	for(i = 0; i < local_readNum; i++){
 		indices[i] = adress_to_write_in_data_by_element[new_rank_phase2[i]];
 		adress_to_write_in_data_by_element[new_rank_phase2[i]] += new_read_size_phase2[i];
 		blocklens[i] = new_read_size_phase2[i];
 		oldtypes[i] = MPI_CHAR;
	 }

	 for(i = 0; i < local_readNum; i++){
	 	assert (indices[i] != (MPI_Aint)NULL);
	 }
	 //Create struct
	 MPI_Type_create_struct(local_readNum, blocklens, indices, oldtypes, &dt_data);
	 MPI_Type_commit(&dt_data);
	 pos=0;
	 res = MPI_Unpack(data_pack, new_data_sz, &pos, MPI_BOTTOM, 1, dt_data, COMM_WORLD);
	 assert(res == MPI_SUCCESS);

	MPI_Type_free(&dt_data);
	MPI_Barrier(COMM_WORLD);

	free(data_pack);
	free(blocklens);
	free(indices);
	free(oldtypes);
	free(new_offset_source_sorted_phase2);
	if (rank == master_job_phase_2)
		fprintf(stderr, "Rank %d :::::[WRITE_ANY_DIM][PACK] Time in packing data for bruck %f seconds\n", rank, MPI_Wtime() - time_count);

	//Free type

	/*
	 * 	In this part we are going to send the data
	 * 	buffer according to the rank it belong
	 * 	and in the sorted order
	 *
	 *  variable we have are:
	 *
	 *  	1) data the buffer = hold the reads
	 *  	2) new_rank_sorted_phase2 vector = hold the rank of the reads
	 *		3) new_offset_dest_sorted_phase2 = hold the offsets of the in the destination file
	 *		4) new_read_size_sorted_phase2 = hold the size of the reads
	 *
	 *
	 *
	 *	The strategy is :
	 *		1) we loop the rank
	 *		2) each rank send to the target rank how much data it's going to send. In order to prepare buffers
	 *		3) for each rank we create datatype of buffered data, read size, and offset
	 *		4) the taget rank recieve the buffered data, the reads size, and the offset vector
	 */

	/****************************
	 * 	BEGIN BRUCK PHASE     	*
	 *****************************/
	//task Beginning of the Bruck phase
	size_t **data_offsets = (size_t **)malloc(sizeof(size_t *)*num_proc);
	int **data_size = (int **)malloc(sizeof(int *)*num_proc);

	time_count = MPI_Wtime();

	bruckWrite(
			rank,
			num_proc,
			local_readNum,
			number_of_reads_by_procs,
			new_rank_phase2,
			buffs_by_procs,
			&data2,
			new_offset_dest_phase2,
			&data_offsets,
			new_read_size_phase2,
			&data_size
		);

	MPI_Barrier(COMM_WORLD);

	free(buffs_by_procs);
	free(new_read_size_phase2);
	free(new_rank_phase2);
	free(new_offset_dest_phase2);

	if (rank == master_job_phase_2)
		fprintf(stderr, "Rank %d :::::[WRITE_ANY_DIM][BRUCK] Time for Bruck reads phases %f seconds\n", rank, MPI_Wtime() - time_count);

	local_readNum = 0;

	for (m = 0; m < num_proc; m++){
		local_readNum += number_of_reads_by_procs[m];
	}

	/*
	 * task: Sort before writing
	 */

	free(new_offset_dest_index_phase2);

	size_t *new_offset_dest_index_phase3 = malloc(sizeof(size_t) * local_readNum);
	//we sort the offset destination
	new_offset_dest_index_phase3[0] = 0;

	for (k = 1; k < local_readNum; k++){
		new_offset_dest_index_phase3[k] = k;
	}

	//task We pack the data before compression
	char  **data_reads_to_sort = malloc(local_readNum * sizeof(char*));
	j=0;
	for(m = 0; m < num_proc; m++)
	{
		int i = 0;
		for(k = 0; k < number_of_reads_by_procs[m]; k++)
		{
			data_reads_to_sort[k + j] = &(data2[m][i]);
			i += data_size[m][k];
			//j += data_size[m][k];
		}
		j += number_of_reads_by_procs[m];
	}

	size_t *data_size_to_sort = malloc(local_readNum * sizeof(size_t));
	j=0;
	for(m = 0; m < num_proc; m++)
	{
		for(k = 0; k < number_of_reads_by_procs[m]; k++)
		{
			data_size_to_sort[k + j] = data_size[m][k];
		}
		free(data_size[m]);
		j += number_of_reads_by_procs[m];
	}
	free(data_size);

	size_t *data_offsets_to_sort = malloc(local_readNum * sizeof(size_t));
	j=0;
	for(m = 0; m < num_proc; m++)
	{
		for(k = 0; k < number_of_reads_by_procs[m]; k++)
		{
			data_offsets_to_sort[k + j] = data_offsets[m][k];
		}
		free(data_offsets[m]);
		j += number_of_reads_by_procs[m];
	}

	if (data_offsets != NULL)
		free(data_offsets);

	free(number_of_reads_by_procs);

	base_arr2 = data_offsets_to_sort;
	qksort(new_offset_dest_index_phase3, local_readNum, sizeof(size_t), 0, local_readNum - 1, compare_size_t);

	size_t* offsets_sorted = malloc(sizeof(size_t)*local_readNum);
	for (k = 0; k < local_readNum; k++){
		offsets_sorted[k] = data_offsets_to_sort[new_offset_dest_index_phase3[k]];
	}

	free(data_offsets_to_sort);

	size_t size_t_buffer_uncompressed = 0;
	for(k = 0; k < local_readNum; k++){
		size_t_buffer_uncompressed += data_size_to_sort[new_offset_dest_index_phase3[k]];
	}

	char *char_buff_uncompressed = malloc(size_t_buffer_uncompressed * sizeof(char) + 1);
	char_buff_uncompressed[size_t_buffer_uncompressed] = 0;
	char *p1 = char_buff_uncompressed;
	size_t q1 = 0;

	for(k = 0; k < local_readNum; k++){
		while ( q1 < data_size_to_sort[new_offset_dest_index_phase3[k]]){
			*p1++ = *data_reads_to_sort[new_offset_dest_index_phase3[k]]++;
			q1++;
		}
		q1=0;
	}


	fprintf (stderr, "[WRITE] rank = %d :: char_buff_uncompressed = %s \n", rank, char_buff_uncompressed);

	free(new_offset_dest_index_phase3);
	free(data_reads_to_sort);

	time_count = MPI_Wtime();

	BGZF *fp;
	fp = calloc(1, sizeof(BGZF));
	int block_length = MAX_BLOCK_SIZE;
	int bytes_written;
	int length = strlen(char_buff_uncompressed);

	fp->open_mode = 'w';
	fp->uncompressed_block_size = MAX_BLOCK_SIZE;
    fp->uncompressed_block = malloc(MAX_BLOCK_SIZE);
    fp->compressed_block_size = MAX_BLOCK_SIZE;
    fp->compressed_block = malloc(MAX_BLOCK_SIZE);
	fp->cache_size = 0;
	fp->cache = kh_init(cache);
	fp->block_address = 0;
	fp->block_offset = 0;
	fp->block_length = 0;
	fp->compress_level = compression_level < 0? Z_DEFAULT_COMPRESSION : compression_level; // Z_DEFAULT_COMPRESSION==-1

	if (fp->compress_level > 9) fp->compress_level = Z_DEFAULT_COMPRESSION;

	const bgzf_byte_t *input = (void *)char_buff_uncompressed;
	int compressed_size = 0;

	if (fp->uncompressed_block == NULL)
	   fp->uncompressed_block = malloc(fp->uncompressed_block_size);

	input = (void *)char_buff_uncompressed;
	block_length = fp->uncompressed_block_size;
	bytes_written = 0;
	uint8_t *compressed_buff =  malloc(strlen(char_buff_uncompressed) * sizeof(uint8_t));

	if (rank == master_job_phase_2)
		fprintf(stderr, "rank %d :::: start loop compression \n", rank);

	while (bytes_written < length) {
		int copy_length = bgzf_min(block_length - fp->block_offset, length - bytes_written);
		bgzf_byte_t* buffer = fp->uncompressed_block;
		memcpy(buffer + fp->block_offset, input, copy_length);
		fp->block_offset += copy_length;
		input += copy_length;
		bytes_written += copy_length;
		//if (fp->block_offset == block_length) {
			//we copy in a temp buffer
			while (fp->block_offset > 0) {
				int block_length;
				block_length = deflate_block(fp, fp->block_offset);

				//is it necessary?
				//if (block_length < 0) break;

				// count = fwrite(fp->compressed_block, 1, block_length, fp->file);
				// we replace the fwrite with a memcopy
				memcpy(compressed_buff + compressed_size, fp->compressed_block, block_length);
	        	compressed_size +=block_length;
	        	fp->block_address += block_length;
			}
	    //}
	}
	if (rank == master_job_phase_2)
		fprintf(stderr, "Rank %d :::::[WRITE_ANY_DIM] Time for compressing %f seconds :::: uncompressed size = %d ::: compression size = %d \n",
			rank, MPI_Wtime() - time_count, length, compressed_size);

	//we compress the neader


	BGZF *fp_header;
	fp_header = calloc(1, sizeof(BGZF));
	uint8_t *compressed_header = NULL;
	int compressed_size_header = 0;

	if (rank == 0) {

		int block_length = MAX_BLOCK_SIZE;
		int bytes_written;
		int length = strlen(header);

		fp_header->open_mode = 'w';
		fp_header->uncompressed_block_size = MAX_BLOCK_SIZE;
		fp_header->uncompressed_block = malloc(MAX_BLOCK_SIZE);
		fp_header->compressed_block_size = MAX_BLOCK_SIZE;
		fp_header->compressed_block = malloc(MAX_BLOCK_SIZE);
		fp_header->cache_size = 0;
		fp_header->block_address = 0;
		fp_header->block_offset = 0;
		fp_header->block_length = 0;
		fp_header->compress_level = compression_level < 0? Z_DEFAULT_COMPRESSION : compression_level; // Z_DEFAULT_COMPRESSION==-1

		if (fp_header->compress_level > 9) fp_header->compress_level = Z_DEFAULT_COMPRESSION;

		const bgzf_byte_t *input = (void *)header;

		if (fp_header->uncompressed_block == NULL)
			fp_header->uncompressed_block = malloc(fp_header->uncompressed_block_size);

		input = (void *)header;
		block_length = fp_header->uncompressed_block_size;
		bytes_written = 0;
		compressed_header =  malloc(strlen(char_buff_uncompressed) * sizeof(uint8_t));

		while (bytes_written < length) {
			int copy_length = bgzf_min(block_length - fp_header->block_offset, length - bytes_written);
			bgzf_byte_t* buffer = fp_header->uncompressed_block;
			memcpy(buffer + fp_header->block_offset, input, copy_length);
			fp_header->block_offset += copy_length;
			input += copy_length;
			bytes_written += copy_length;
			//if (fp->block_offset == block_length) {
			//we copy in a temp buffer
			while (fp_header->block_offset > 0) {
				int block_length;
				block_length = deflate_block(fp_header, fp_header->block_offset);

				//is it necessary?
				//if (block_length < 0) break;

				// count = fwrite(fp->compressed_block, 1, block_length, fp->file);
				// we replace the fwrite with a memcopy
				memcpy(compressed_header + compressed_size_header, fp_header->compressed_block, block_length);
				compressed_size_header +=block_length;
				fp_header->block_address += block_length;
			}
			//}
		}
	}
	
	kh_destroy(cache, fp->cache);
	
	free(char_buff_uncompressed);
	size_t compSize = compressed_size;

	/*
	 * We write results of compression
	 */
	MPI_Barrier(COMM_WORLD);
	size_t write_offset = 0;

	MPI_Offset * y  = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
	MPI_Offset * y2 = (MPI_Offset *) calloc(num_proc+1, sizeof(MPI_Offset));

	MPI_Gather(&compSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);

	//now we make a cumulative sum
	int i1 = 0;

	if (rank ==0){
		for (i1 = 1; i1 < (num_proc + 1); i1++) {
			y2[i1] = y[i1-1];
			}

		for (i1 = 1; i1 < (num_proc +1); i1++) {
			y2[i1] = y2[i1-1] + y2[i1];
		}

		for (i1 = 0; i1 < (num_proc +1); i1++) {
			y2[i1] = y2[i1] + write_offset + compressed_size_header;
		}

	}

	MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);
	// we create the path where to write for collective write
	path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
	sprintf(path, "%s/%s.gz", output_dir, chrName);

	if(!rank)
		fprintf(stderr, "Rank %d :::::[WRITE_ANY_DIM] Opening the file %s \n", rank, path );


	//task FINE TUNING FINFO FOR WRITING OPERATIONS
	/*
	MPI_Info_set(finfo,"striping_factor","128");
	MPI_Info_set(finfo,"striping_unit","1610612736"); //1G striping
	MPI_Info_set(finfo,"nb_proc","128");
	MPI_Info_set(finfo,"cb_nodes","128");
	MPI_Info_set(finfo,"cb_block_size","1610612736"); // 4194304 MBytes - should match FS block size
	MPI_Info_set(finfo,"cb_buffer_size","1610612736"); // 128 MBytes (Optional)
	*/

	ierr = MPI_File_open(COMM_WORLD, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);

	if (ierr) {
		fprintf(stderr, "Rank %d :::::[WRITE_ANY_DIM] failed to open %s.\nAborting.\n\n", rank, path);
		MPI_Abort(COMM_WORLD, ierr);
		exit(2);
	}
	else{
		if(!rank)fprintf(stderr, "Rank %d :::[WRITE_ANY_DIM] %s.bam successfully opened\n", rank, chrName);
	}

	time_count = MPI_Wtime();

	if (rank == master_job_phase_2 ) {
		fprintf(stderr, "Rank %d :::::[WRITE_ANY_DIM] we write the header \n", rank);
		MPI_File_write(out, compressed_header, compressed_size_header, MPI_BYTE, MPI_STATUS_IGNORE);
	}
	free(compressed_header);

	MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
	MPI_File_write_all(out, compressed_buff, (size_t)compSize, MPI_BYTE, &status);

	//task FINE TUNING FINFO BACK TO READING OPERATIONS
	/*
	MPI_Info_set(finfo,"striping_factor","128");
	MPI_Info_set(finfo,"striping_unit","2684354560"); //1G striping
	MPI_Info_set(finfo,"nb_proc","128");
	MPI_Info_set(finfo,"cb_nodes","128");
	MPI_Info_set(finfo,"cb_block_size","2684354560"); // 4194304 MBytes - should match FS block size
	MPI_Info_set(finfo,"cb_buffer_size","2684354560"); // 128 MBytes (Optional)
	*/

	free(compressed_buff);

	if (rank == master_job_phase_2)
		fprintf(stderr, "Rank %d :::::[WRITE_ANY_DIM] Time for chromosome %s writing %f seconds\n\n\n",
				rank, chrName, MPI_Wtime()-time_count);

	free(fp->uncompressed_block);
	free(fp->compressed_block);
	free_cache(fp);
	free(fp);

	if (rank == 0){
		free(fp_header->uncompressed_block);
		free(fp_header->compressed_block);
		free_cache(fp_header);
	}
	free(fp_header);


	MPI_File_close(&out);
	free(path);
	for(m = 0; m < num_proc; m++)
	{
		if (data2[m]) free(data2[m]);
	}
	if (data2) free(data2);

	free(offsets_sorted);
	free(new_read_size_sorted_phase3);
	free(data_size_to_sort);
	free(y);
	free(y2);

}





