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
#include "bgzf2.h"
#include "bgzf2.c"
#include "mpiSortUtils.h"
#include "parallelBitonicSort2.h"
#include "mergeSort.h"
#include "libgen.h"

#ifdef HAVE_HTSLIB
#include "sam.h"
#include "hfile.h"
#include "hfile.c"
#include "hts.h"
#endif

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

void bruckWrite4(
                int rank,
                int num_proc,
                size_t local_readNum,
                size_t* number_of_reads_by_procs,
                int *new_rank,
                int *new_reads_size,
                int ***read_size
                ){

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
		size_t previous_local_readNum,
		int uniq_chr,
		int write_format,
		int merge,
		char file_name_sorted[]
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
		//The data in which what we read will be kept during Bruck
		//char **data2 = malloc( total_num_proc * sizeof(char*));
		
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
			new_pbs_offset_dest_phase1 	= calloc( (local_readNum - tmp2), sizeof(size_t));
			new_pbs_read_size_phase1 	= calloc( (local_readNum - tmp2), sizeof(int));
			new_pbs_dest_rank_phase1 	= calloc( (local_readNum - tmp2), sizeof(int));
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
			 */
			for (j = 0; j < (local_readNum - tmp2); j++){
				assert(new_pbs_read_size_phase1[j]     != 0);
				assert(new_pbs_offset_source_phase1[j] != 0);
				assert(new_pbs_offset_dest_phase1[j]   != 0);
				assert(new_pbs_dest_rank_phase1[j]     < dimensions);
				assert(new_pbs_orig_rank_off_phase1[j] < dimensions);
			}
			
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

		int *new_local_reads_sizes_sorted_bruck 	= malloc(previous_local_readNum * sizeof(int));
		int *new_local_reads_dest_rank_sorted_bruck   	= malloc(previous_local_readNum * sizeof(int));
		size_t *new_local_offset_destination_bruck 	= malloc(previous_local_readNum * sizeof(size_t));
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


		size_t **dest_offsets 		= malloc(sizeof(size_t *) * dimensions);
		size_t **local_source_offsets 	= malloc(sizeof(size_t *) * dimensions);

		int **read_size 		= malloc(sizeof(int *) * dimensions);
		int **dest_rank 		= malloc(sizeof(int *) * dimensions);


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
				new_local_reads_sizes_sorted_bruck[k + j] 	= read_size[m][k];
				new_local_offset_destination_bruck[k + j] 	= dest_offsets[m][k];
			}
			free(local_source_offsets[m]);
			free(dest_rank[m]);
			free(read_size[m]);
			free(dest_offsets[m]);
			j += number_of_reads_by_procs[m];

		}
		assert( j == previous_local_readNum );

		for(k = 1; k < previous_local_readNum; k++){
			assert(new_local_offset_source_sorted_bruck[k-1] < new_local_offset_source_sorted_bruck[k]);
		}


		time_count = MPI_Wtime();
		//init indices for qksort
		size_t *coord_index = malloc(previous_local_readNum * sizeof(size_t));

		for(j = 0; j < previous_local_readNum; j++){
			coord_index[j] = j;
		}

		base_arr2 = new_local_offset_source_sorted_bruck;

		// Replace qksort with merge sort for stability
		//qksort(coord_index, previous_local_readNum, sizeof(size_t), 0, previous_local_readNum - 1, compare_size_t);
		MergeSortMain( coord_index, previous_local_readNum);

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
		 */
		

		for(j = 0; j < previous_local_readNum - 1; j++){

			assert(new_local_offset_source_sorted_bruck2[j] < new_local_offset_source_sorted_bruck2[j+1]);
		}
		

		free(new_local_offset_source_sorted_bruck);
		free(new_local_reads_sizes_sorted_bruck);
		free(new_local_offset_destination_bruck);
		free(new_local_reads_dest_rank_sorted_bruck);
		free(coord_index);

		//malloc_trim(0);

		if (rank == master_job_phase_2)
			fprintf(stderr, "Rank %d :::::[WRITE][LOCAL SORT] Time =  %f seconds\n", rank, MPI_Wtime() - time_count);

		/*
		 *
		 * FOR DEBUG
		 *
		 */
		
		for(k = 1; k < previous_local_readNum; k++)
		{
			assert((new_local_offset_source_sorted_bruck2[k] - start_offset_in_file) <= strlen(data));
			assert(new_local_offset_source_sorted_bruck2[k-1]    < new_local_offset_source_sorted_bruck2[k]);
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
		

		/*  Before moving the data with bruck
 		*   We are going to compute the buffer 
 		*   size needed for the read in order to allocate exactly 
 		*   data2
 		*
 		*  we use 
 		*  new_local_reads_dest_rank_sorted_bruck2
 		*  new_local_reads_sizes_sorted_bruck2
 		*/

		time_count = MPI_Wtime();
		size_t *sz_pack_by_proc = calloc(num_proc,sizeof(size_t));
                size_t *number_of_reads_by_procs_tmp = calloc( dimensions, sizeof(size_t));
                for(m = 0; m < previous_local_readNum; m++) 
			number_of_reads_by_procs_tmp[new_local_reads_dest_rank_sorted_bruck2[m]]++;
                count6 = 0;
                for(m = 0; m < dimensions; m++){
                                count6 += number_of_reads_by_procs_tmp[m];
                }
		assert( count6 == previous_local_readNum);

                int **read_size_2                 = malloc(sizeof(int *) * dimensions);
                bruckWrite4(
                                rank,
                                dimensions,
                                count6,
                                number_of_reads_by_procs_tmp,
                                new_local_reads_dest_rank_sorted_bruck2,
                                new_local_reads_sizes_sorted_bruck2,
                                &read_size_2
                );
			
		count5=0;
                j=0;
                for(m = 0; m < num_proc; m++){
                        for(k = 0; k < number_of_reads_by_procs_tmp[m]; k++)  
				sz_pack_by_proc[m] += read_size_2[m][k];
                        free(read_size_2[m]);
                        j += number_of_reads_by_procs_tmp[m];

                }
                assert( j == previous_local_readNum );
		free(number_of_reads_by_procs_tmp);
		free(read_size_2);	
		
		if (rank == master_job_phase_2)
                        fprintf(stderr, "Rank %d :::::[WRITE][COMPUTE BUFFER SIZE] Time =  %f seconds\n", rank, MPI_Wtime() - time_count);


		// Now we are ready to exchange the data	
		// We have the vectors
		// new_local_reads_sizes_sorted_bruck2
		// new_local_reads_dest_rank_sorted_bruck2[k]
		// new_local_offset_destination_bruck2[k] 
		// new_local_offset_source_sorted_bruck2[k]
		
		//we compute if we need to use bruck several times or not
		//we set the limits to 2GB
		//we gonna build pack of 2 GB size and do the Bruck in several times
 		//we compute number of packs we have to send 

		//#define STRIPING_UNIT "1610612736"  // 1.5GB
		//#define STRIPING_UNIT "2147483648"  // 2GB

		//size_t pack_size = 2147483648;
		size_t pack_size_tmp = 1024*1024*1024; //1GB
                int chunk_data_number = 0;
                size_t tmp_2 = 0;

                 for (k = 0; k < previous_local_readNum; k++){
                        assert (new_local_reads_sizes_sorted_bruck2[k] != 0);
                        tmp_2 += new_local_reads_sizes_sorted_bruck2[k];
                }

                int pack_number = (tmp_2 / pack_size_tmp) + 1;
                int max_number_of_packs = 0;
                //due to load balancing all ranks must have the same number of pack 
                MPI_Allreduce(&pack_number, &max_number_of_packs, 1, MPI_LONG_LONG_INT, MPI_MAX, split_comm);
                pack_number = max_number_of_packs;
                //according to the max number we compute a new pack_size per rank
                size_t pack_size = 0;
                pack_size = tmp_2 / max_number_of_packs;
                size_t **data_offsets   = malloc(sizeof(size_t *) * num_proc);
                int **data_size         = malloc(sizeof(int *) * num_proc);
                if (rank == master_job_phase_2)
                	fprintf (stderr, "Rank %d :::: [WRITE][COMPUTE PACKs] Number of packs buffer = %d \n", rank, pack_number);
                //
		//hold the index of each reads in data_reads_to_sort 
		char **data_reads_to_sort      = calloc(previous_local_readNum, sizeof(char *));		
		char **data2 = malloc(( pack_number * num_proc) * sizeof(char *));
	
		int *data_size_to_sort          	= calloc(previous_local_readNum, sizeof(int));
		assert(data_size_to_sort);
                size_t *data_offsets_to_sort    	= calloc(previous_local_readNum, sizeof(size_t));
		assert(data_offsets_to_sort);
		size_t *number_of_reads_by_procs2       = calloc( dimensions , sizeof(size_t));
		assert(number_of_reads_by_procs2);
		
		/*
		for ( k = 0; k < num_proc ; k++){
			
			 data2[k] = malloc((sz_pack_by_proc[k] + 1)*sizeof(char));			
			 data2[k][sz_pack_by_proc[m]] = 0;
			 number_of_reads_by_procs[k] = 0; 
		}
		*/
	
		free(sz_pack_by_proc);
		//In a for loop we send each number of pack reads
		int pack_index = 0;
		int packs = 0;
		char *p = NULL;
		char *q = NULL;
		size_t sz_orig_data = strlen(data);
		size_t k1 = 0;
		tmp_2 = 0;
		size_t start_index = 0;
		size_t last_index = 0;
		size_t num_read_to_pack = 0;
		size_t size_to_pack = 0;
		size_t number_of_reads_recieved = 0;
		size_t total_size_recv = 0;
                int *new_local_reads_sizes_sorted_bruck3 = NULL;
                int* new_local_reads_dest_rank_sorted_bruck3 = NULL;
                size_t *new_local_offset_destination_bruck3 = NULL;
                size_t *new_local_offset_source_sorted_bruck3 = NULL;
		size_t *previous_sz_pack_by_proc = calloc(num_proc,sizeof(size_t));
		double total_time_for_bruck = MPI_Wtime();
	
		for ( k1 = 0; k1 < previous_local_readNum; k1++){
 			// we build new vector which point to
 			// previous vector at the right start
 			//size of the new_vector are number_reads_in_pack 
 			tmp_2 += new_local_reads_sizes_sorted_bruck2[k1];

			if (tmp_2 < pack_size && (k1 < (previous_local_readNum -1))) continue;		
			else{
                 		
				if ( k1 == (previous_local_readNum - 1)){
					size_to_pack = tmp_2;
                                        last_index = k1;
                                	num_read_to_pack =  last_index - start_index + 1;
                                 }
                                 else{
					size_to_pack = tmp_2 - new_local_reads_sizes_sorted_bruck2[k1];
                                        last_index = k1 - 1;
					num_read_to_pack =  last_index - start_index + 1;
                                 }
	
				size_t current_pack_readNum = num_read_to_pack;
				/*
 				*
 				* FOR DEBUG
 				*
                               	if (rank == 0){
	 	
					 fprintf(stderr, "rank %d : pack_index = %d : chunk:  start = %zu ::: end =%zu \n", rank, pack_index, start_index, last_index);	
					 fprintf(stderr, "rank %d : pack_index = %d : previous_local_readNum = %zu \n", rank, pack_index, previous_local_readNum);
                                         fprintf(stderr, "rank %d : pack_index = %d : current pack readnum = %zu \n", rank, pack_index, current_pack_readNum);
					 fprintf(stderr, "rank %d : pack_index = %d : sz_orig_data = %zu \n", rank, pack_index, sz_orig_data);	
					 fprintf(stderr, "rank %d : pack_index = %d : size_to_pack = %zu \n", rank, pack_index, size_to_pack);

				}
				*/	
				new_local_reads_sizes_sorted_bruck3 	= calloc(num_read_to_pack, sizeof(int));
				new_local_reads_dest_rank_sorted_bruck3 = calloc(num_read_to_pack, sizeof(int));
				new_local_offset_destination_bruck3 	= calloc(num_read_to_pack, sizeof(size_t));
				new_local_offset_source_sorted_bruck3 	= calloc(num_read_to_pack, sizeof(size_t));	
				
				for ( k = 0; k < num_read_to_pack; k++){
					new_local_reads_sizes_sorted_bruck3[k] 		= new_local_reads_sizes_sorted_bruck2[start_index + k]; 
					new_local_reads_dest_rank_sorted_bruck3[k] 	= new_local_reads_dest_rank_sorted_bruck2[start_index + k];
					new_local_offset_destination_bruck3[k]		= new_local_offset_destination_bruck2[start_index + k];
					new_local_offset_source_sorted_bruck3[k]        = new_local_offset_source_sorted_bruck2[start_index + k];
				
					assert( new_local_offset_source_sorted_bruck3[k] > 0);
					assert( new_local_offset_destination_bruck3[k] != 0);
					assert( new_local_reads_sizes_sorted_bruck3[k] != 0);
					assert ((new_local_offset_source_sorted_bruck3[k] - start_offset_in_file) < sz_orig_data);
				}
				
				start_index = k1;
				/*
				*new_local_reads_sizes_sorted_bruck3[0] = (new_local_reads_sizes_sorted_bruck2 + start_index);
				*new_local_reads_dest_rank_sorted_bruck3[0] = (new_local_reads_dest_rank_sorted_bruck2 + start_index);
				*new_local_offset_destination_bruck3[0] = (new_local_offset_destination_bruck2 + start_index);
				*new_local_offset_source_sorted_bruck3[0] = (new_local_offset_source_sorted_bruck2 + start_index);
				*/
				//we pack 
				size_t i;
				size_t new_data_sz = 0;
				char *data_pack;
				char *tmp_tab;

				//we compute the size of data_pack
				for (k = 0; k < current_pack_readNum; k++){
					new_data_sz += new_local_reads_sizes_sorted_bruck3[k];
				}
			
				data_pack = malloc(new_data_sz +1);
				data_pack[new_data_sz] = 0;

				q = data;
				p = data_pack;
				size_t offset_in_data = 0;
				int pos = 0;
				
				//we copy elements from data in data_pack
				for (k=0; k < (current_pack_readNum); k++){
					pos = 0;
					offset_in_data = new_local_offset_source_sorted_bruck3[k] - start_offset_in_file;
					assert(offset_in_data < sz_orig_data);
					q = data + offset_in_data;
					assert( *q != 0);
					while (*q && (pos < new_local_reads_sizes_sorted_bruck3[k])) { assert(*q != 0); *p=*q; q++; p++; pos++; }
				}

				/*
				if (rank == 0)
                                       	fprintf(stderr, "rank %d ::: pack_index = %d ::: strlen(data_pack) = %zu \n", rank, pack_index, strlen(data_pack));
				*/

				assert( data_pack[0] != 0);
				assert( strlen(data_pack)> 0);
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

				for(m = 0; m < num_proc; m++) number_of_reads_by_procs2[m] = 0;

				for(m = 0; m < current_pack_readNum; m++)
				{
					buffs_by_procs2[new_local_reads_dest_rank_sorted_bruck3[m]] += new_local_reads_sizes_sorted_bruck3[m];
					number_of_reads_by_procs2[new_local_reads_dest_rank_sorted_bruck3[m]]++;
				}

				for(m = 0; m < num_proc; m++)
				{
					buffs_by_procs[(rank + m)%num_proc] = buffs_by_procs2[(rank - m + num_proc)%num_proc];
				}

				free(buffs_by_procs2);

				//Allocate data and initialization
				char **data3                            = malloc( num_proc * sizeof(char*));
				for(m = 0; m < num_proc; m++){
					data3[m] = (char*)malloc(buffs_by_procs[m]*sizeof(char) +1);
					assert(data3[m]);
					data3[m][buffs_by_procs[m]] = 0;
				}

				//Variable for datatype struct
				MPI_Aint *indices 		= malloc(current_pack_readNum * sizeof(MPI_Aint));
				int *blocklens    		= malloc(current_pack_readNum * sizeof(int));
				MPI_Datatype *oldtypes 		= malloc(current_pack_readNum * sizeof(MPI_Datatype));

				MPI_Aint adress_to_write_in_data_by_element[num_proc];
				for(i = 0; i < num_proc; i++){
					MPI_Get_address(data3[(rank-i+num_proc)%num_proc], &adress_to_write_in_data_by_element[(rank+i)%num_proc]);
				}

				for(i = 0; i < current_pack_readNum; i++){
					indices[i] = adress_to_write_in_data_by_element[new_local_reads_dest_rank_sorted_bruck3[i]];
					assert (indices[i] != (MPI_Aint)NULL);
					adress_to_write_in_data_by_element[new_local_reads_dest_rank_sorted_bruck3[i]] += new_local_reads_sizes_sorted_bruck3[i];
					blocklens[i] = new_local_reads_sizes_sorted_bruck3[i];
					oldtypes[i] = MPI_CHAR;
			 	}

				//Create struct
				MPI_Type_create_struct(current_pack_readNum, blocklens, indices, oldtypes, &dt_data);
				MPI_Type_commit(&dt_data);
				pos=0;
				res = MPI_Unpack(data_pack, new_data_sz, &pos, MPI_BOTTOM, 1, dt_data, COMM_WORLD);
				assert(res == MPI_SUCCESS);
				MPI_Type_free(&dt_data);

				free(data_pack);
				free(blocklens);
				free(indices);
				free(oldtypes);

				if (rank == master_job_phase_2)
					fprintf(stderr, "Rank %d :::::[WRITE][DATA PACK %d] : Time =  %f seconds\n", rank, packs, 
							MPI_Wtime() - time_count);

				/****************************
				 * 	BEGIN BRUCK PHASE     	*
				 *****************************/

				size_t **data_offsets2 		= malloc(sizeof(size_t *) * num_proc);
				int **data_size2 	   	= malloc(sizeof(int *) * num_proc);

				time_count = MPI_Wtime();

				count6 = 0;
				for(m = 0; m < dimensions; m++){
					count6 += number_of_reads_by_procs2[m];
				}
					
				assert( count6 == current_pack_readNum );

				bruckWrite(
					rank,
					num_proc,
					current_pack_readNum,
					number_of_reads_by_procs2,
					new_local_reads_dest_rank_sorted_bruck3,
					buffs_by_procs,
					&data3,
					new_local_offset_destination_bruck3,
					&data_offsets2,
					new_local_reads_sizes_sorted_bruck3,
					&data_size2
				);
				if (rank == master_job_phase_1)
					fprintf(stderr, "Rank %d :::::[WRITE][BRUCK PACK %d] :: Time   = %f s \n", rank, packs, MPI_Wtime() - time_count);
	
				for(m = 0; m < num_proc; m++)
				{
					size_t buff_test  = strlen(data3[m]);
					size_t buff_test2 = 0;
					size_t u = 0;
					for (u = 0; u < number_of_reads_by_procs2[m]; u++){
						buff_test2 += data_size2[m][u];
						assert(data_offsets2[m][u] != 0);
						assert(data_size2[m][u] != 0);
					}
					assert (buff_test2 == buff_test);
				}

				free(buffs_by_procs);
				/*
			 	* GET DATA AFTER BRUCK
			 	*
			 	*/

				// data3 will hold the result of all the bruck
				// we copy data3 in data2
				//for assert
				j = 0;
				//we copy contains of data3 in data2
				size_t sz_data3 = 0;
				size_t sz_data2 = 0;
				for(m = 0; m < num_proc; m++){
					
					data2[pack_index] = malloc( (strlen(data3[m]) + 1)*sizeof(char));
                         		data2[pack_index][strlen(data3[m])] = 0;
					// fist we update data2
					sz_data3 = strlen(data3[m]);
					char *u = data2[pack_index];
					char *v = data3[m];
					while(*v){*u = *v; u++; v++;}
					total_size_recv += sz_data3;
					free(data3[m]);
					
					
					//then we update
					// data_reads_to_sort
					// data_size_to_sort
					// data_offsets_to_sort
					/*
 					 * FOR DEBUG
 					 *
						for(k = 0; k < number_of_reads_by_procs[m]; k++)
                                        	{
							assert(data_size_to_sort[k] != 0 );
							assert(data_offsets_to_sort[k] != 0);
						}
					}
					*/
					i = 0;
					//**data_reads_to_sort = malloc ( number_of_reads_by_procs2[m] * sizeof(char*));

					for(k = 0; k < number_of_reads_by_procs2[m]; k++)
					{
						// we start to fill up the vector after the previous 
						// number of reads in all previous buffer
						// memcopy should be faster
						
						data_size_to_sort[number_of_reads_recieved + k] = data_size2[m][k];
                                                assert( data_size_to_sort[number_of_reads_recieved + k] != 0 );

                                                data_offsets_to_sort[number_of_reads_recieved + k] = data_offsets2[m][k];
                                                assert(data_offsets_to_sort[number_of_reads_recieved + k ] != 0);
						
						data_reads_to_sort[number_of_reads_recieved + k] = &(data2[pack_index][i]);		
						i += data_size2[m][k];
					}
						
					previous_sz_pack_by_proc[m] = sz_data3;		
					j +=  number_of_reads_by_procs2[m];
					number_of_reads_recieved += number_of_reads_by_procs2[m];
					free(data_size2[m]);
                                       	free(data_offsets2[m]);	
					number_of_reads_by_procs[m] += number_of_reads_by_procs2[m];
					//we reset number_of_reads_by_procs2[m]
					number_of_reads_by_procs2[m] = 0;
					pack_index++;
				}
				
				//assert(j == current_pack_readNum);
				/*
 				 * FOR DEBUG
 				 *
 				 *
				if (rank == 0){
					fprintf(stderr,  "Rank %d ::: pack :: %d :: total reads recieved =  %zu \n", rank, pack_index, number_of_reads_recieved);
					fprintf(stderr,  "Rank %d ::: pack :: %d :: total size recieved =  %zu \n", rank, pack_index, total_size_recv);
				}
				*/
				free(data3);
				
				free(new_local_reads_sizes_sorted_bruck3);
                               	free(new_local_reads_dest_rank_sorted_bruck3);
                               	free(new_local_offset_destination_bruck3);
				free(new_local_offset_source_sorted_bruck3);
				
                         	free(data_size2);
				free(data_offsets2);
				tmp_2=0;
				
			}//end if tmp_2 > pack_size	

			packs++;		
		} //end for loop on previous_local_num_read
	
		free(previous_sz_pack_by_proc);
		assert(number_of_reads_recieved == previous_local_readNum);
		/*
 		 *	FOR DEBUG
 		 *
 		 *
		size_t total_reads = 0;
		MPI_Allreduce(&previous_local_readNum, &total_reads, 1, MPI_LONG_LONG_INT, MPI_SUM, split_comm);
		fprintf(stderr,  "Rank %d :: total reads recieved =  %zu \n", rank, total_reads);

		if (rank == 0){
	                fprintf(stderr,  "Rank %d :: total reads recieved =  %zu \n", rank, number_of_reads_recieved);
			fprintf(stderr,  "Rank %d :: previous local read num =  %zu \n", rank, previous_local_readNum );
			fprintf(stderr,  "Rank %d :: size of data_reads_to_sort =  %zu \n", rank, total_size_recv );
		}
		assert(number_of_reads_recieved == previous_local_readNum);
		size_t tmp5 = 0;
		for(k = 0; k < num_proc; k++)
                 {  
			tmp5 += number_of_reads_by_procs[k];
		}
                   
		if (rank == 0){
                	fprintf(stderr,  "Rank %d :: tmp5 =  %zu \n", rank, tmp5);
			assert(tmp5 == number_of_reads_recieved);

		}
		*/
		 if (rank == 0)
                        fprintf(stderr,  "Rank %d :::::[WRITE][BRUCK END] finish bruck total time %f s\n", rank, MPI_Wtime() - total_time_for_bruck );

		//now free
		//free(data3);
		free(number_of_reads_by_procs2);	
		if (uniq_chr) free(data);

		free(new_local_reads_sizes_sorted_bruck2);
		free(new_local_offset_destination_bruck2);
		free(new_local_reads_dest_rank_sorted_bruck2);
		/*
		for ( k = 0; k < previous_local_readNum; k++) 
			assert(*data_reads_to_sort[k] == 'D');	
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
		
		time_count=MPI_Wtime();
		
		size_t *new_offset_dest_index_phase3	= malloc(sizeof(size_t) * previous_local_readNum);
		for (k = 0; k < previous_local_readNum; k++){
			new_offset_dest_index_phase3[k] = k;
		}
		//task SORT OFFSET DESTINATION
		//now we sort new_offset_dest_phase2


		base_arr2 = data_offsets_to_sort;

		//replace qksort with mergesort
		//qksort(new_offset_dest_index_phase3, previous_local_readNum, sizeof(size_t), 0, previous_local_readNum - 1, compare_size_t);
		MergeSortMain( new_offset_dest_index_phase3, previous_local_readNum );		


		free(data_offsets_to_sort);


		if (rank == master_job_phase_2)
                        fprintf(stderr, "Rank %d :::::[WRITE][BEFORE COMPRESSION] local sort in %f seconds \n",
                                        rank, MPI_Wtime()-time_count);
                        
		time_count = MPI_Wtime();

		size_t size_t_buffer_uncompressed = 0;
		for(k = 0; k < previous_local_readNum; k++){
			size_t_buffer_uncompressed += data_size_to_sort[new_offset_dest_index_phase3[k]];
		}

		assert( size_t_buffer_uncompressed == total_size_recv);	

		char *char_buff_uncompressed = calloc(size_t_buffer_uncompressed + 1, sizeof(char));
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
		

		assert(strlen(char_buff_uncompressed) == total_size_recv);
		free(new_offset_dest_index_phase3);

		free(data_reads_to_sort);
		for ( k = 0; k < (num_proc * pack_number) ; k++) free(data2[k]);
                free(data2);
		free(data_size_to_sort);
	
		if (rank == master_job_phase_2)
                	fprintf(stderr, "Rank %d :::::[WRITE][BEFORE WRITING] internal copy in %f seconds \n",
                                        rank, MPI_Wtime()-time_count);



		if (write_format == 2){
			
			/*
                         * We write result in SAM format
                         */
			time_count = MPI_Wtime();
			int file_exist = 1; //incase of merge test if the file exist
                        size_t write_offset = 0;
			size_t samSize = strlen(char_buff_uncompressed);				
			size_t size_header = strlen(header);
			size_t tmp_size_buffer = 1024*1024*1024;
			size_t tmp_size_buffer2 = tmp_size_buffer;
			
                        MPI_Offset * y  = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
                        MPI_Offset * y2 = (MPI_Offset *) calloc(num_proc+1, sizeof(MPI_Offset));

                        MPI_Gather(&samSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);

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
                                     y2[i1] = y2[i1] + write_offset + size_header;
                              }
                        }
                        MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);
                                 
			
			if (!merge){
				path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
                        	sprintf(path, "%s/%s.sam", output_dir, chrName);
				/* we write the header */
                        	ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);
			}
                        else{
				path = (char*)malloc((strlen(output_dir) + strlen(file_name_sorted) + 40) * sizeof(char));
					
				sprintf(path, "%s/%s", output_dir, file_name_sorted);

				//we test if the file exist
				file_exist = access(path, F_OK);
				ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY|MPI_MODE_CREATE|MPI_MODE_APPEND, finfo, &out);
					                        	
			}

                        if (ierr) {
                                fprintf(stderr, "Rank %d failed to open %s.\nAborting.\n\n", rank, path);
                                MPI_Abort(COMM_WORLD, ierr);
                                exit(2);
                        }
                        else{
                                if(!rank && merge == 0)fprintf(stderr, "Rank %d :::::[WRITE][OPEN SAM RESULTS] %s.sam successfully opened\n", rank, chrName);
                        	if(!rank && merge == 1)fprintf(stderr, "Rank %d :::::[WRITE][OPEN SAM RESULTS] %s_sorted.sam successfully opened\n", rank, file_name);
			}

                        time_count = MPI_Wtime();

			MPI_Offset file_size;
			MPI_File_get_size(out, &file_size); 

                        if (rank == master_job_phase_2 && (merge == 0)) {
                                MPI_File_write(out, header, size_header, MPI_CHAR, MPI_STATUS_IGNORE);
                        }

			if (rank == master_job_phase_2 && (merge == 1) && (file_size == 0)) {
                                MPI_File_write(out, header, size_header, MPI_CHAR, MPI_STATUS_IGNORE);
                        }

			//case of merge the file already exist
                        if (merge && file_exist == 0) write_offset += (file_size - size_header);
			//case of merge the file don't exist
			if  (merge && file_exist > 0) write_offset += file_size;

			MPI_Barrier(COMM_WORLD);
			/* we write the SAM */
			if (samSize < tmp_size_buffer){
                                MPI_File_set_view(out, write_offset, MPI_CHAR, MPI_CHAR, "native", finfo);
                                MPI_File_write_all(out, char_buff_uncompressed, samSize, MPI_CHAR, &status);
                        }
                        //we write by block of 1gb
                        else {
                        	char *buff_tmp = char_buff_uncompressed;
                                size_t tmp10 = 0;
                                int block_tmp = 0;
                                int count_status = 0;
                                int error_status = 0;
                                while (tmp_size_buffer2 > 0){
                                	MPI_File_set_view(out, write_offset, MPI_CHAR, MPI_CHAR, "native", finfo);
                                	MPI_File_write_all(out, buff_tmp, tmp_size_buffer2, MPI_CHAR, &status);
                                        MPI_Get_count(&status, MPI_CHAR, &count_status);
                                        assert(count_status == tmp_size_buffer2);
                                        buff_tmp += tmp_size_buffer2;
                                        tmp10 += tmp_size_buffer2;
                                        write_offset += tmp_size_buffer2;
                                        block_tmp++;
                                        if ( (samSize - tmp10) > tmp_size_buffer ) tmp_size_buffer2 = tmp_size_buffer;
                                        else tmp_size_buffer2 = (samSize - tmp10);
                                }
                        }
			
			MPI_File_close(&out);
                        
                        if (rank == master_job_phase_2)
                        	fprintf(stderr, "Rank %d :::::[WRITE][WRITING SAM] Time for chromosome %s writing %f seconds\n",
                        	        rank, chrName, MPI_Wtime()-time_count);
                        
			free(char_buff_uncompressed);
	        	free(y);
	                free(y2);
                        MPI_Barrier(COMM_WORLD);
		}
                /*
                 *  BGZF Case
                 *
                 */
                if (write_format == 0){
      

			int file_exist = 1;
                        if (!merge){
	                        path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
                                sprintf(path, "%s/%s.gz", output_dir, chrName);
                                ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);
                        }
                        else{
                                path = (char*)malloc((strlen(output_dir) + strlen(file_name_sorted) + 40) * sizeof(char));
                                sprintf(path, "%s/%s", output_dir, file_name_sorted);

                                file_exist = access(path, F_OK);
                                ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY|MPI_MODE_CREATE|MPI_MODE_APPEND, finfo, &out);
                        }

                        if (ierr) {
                                fprintf(stderr, "Rank %d :::::[WRITE][BGZF RESULTS] failed to open %s.\nAborting.\n\n", rank, path);
                                MPI_Abort(COMM_WORLD, ierr);
                                exit(2);
                        }
                        else{
                                if((rank == master_job_phase_2) && merge == 0)fprintf(stderr, "Rank %d :::::[WRITE][BGZF RESULTS] %s.gz successfully opened\n", rank, chrName);
                                if((rank == master_job_phase_2) && merge == 1)fprintf(stderr, "Rank %d :::::[WRITE][BGZF RESULTS] %s_sorted.gz successfully opened\n", rank, file_name);
                        }
	
			if (rank == master_job_phase_2)
                        	fprintf(stderr, "Rank %d :::::[WRITE][BEFORE COMPRESSION] start compression \n", rank);
			
			
			/*
		 	* COMPRESSION PART
		 	*
		 	*/
		
			time_count = MPI_Wtime();
			//the buffer we will use to compress
			size_t tmp_size_buffer = 1024*1024*1024;
			//char *tmp_buffer =malloc(tmp_size_buffer);
			//tmp_buffer[tmp_size_buffer]=0;
			char *p6 = char_buff_uncompressed;
			//char *p2 = tmp_buffer;
			uint8_t *compressed_buff =  malloc((strlen(char_buff_uncompressed))* sizeof(uint8_t));
			assert(compressed_buff);	
			int m6 = 0;
			size_t compressed_size =0;
			BGZF2 *fp;
                	fp = calloc(1, sizeof(BGZF2));	
			while (*p6){

				char *tmp_buffer =malloc(tmp_size_buffer + 1);
				tmp_buffer[tmp_size_buffer]=0;
				char *p7 = tmp_buffer;
				size_t counter_tmp = 0;
				while (*p6 && counter_tmp < tmp_size_buffer) {*p7=*p6; p6++;p7++; counter_tmp++;} 

				time_count = MPI_Wtime();

				//BGZF *fp;
				//fp = calloc(1, sizeof(BGZF));
				int block_length = MAX_BLOCK_SIZE;
				int bytes_written;
				//size_t length = strlen(char_buff_uncompressed);
				int length = counter_tmp;

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

				//const bgzf_byte_t *input = (void *)char_buff_uncompressed;
				//int compressed_size = 0;

				if (fp->uncompressed_block == NULL)
		   			fp->uncompressed_block = malloc(fp->uncompressed_block_size);

				//input = (void *)char_buff_uncompressed;
				const bgzf_byte_t *input = (void *)tmp_buffer;
				block_length = fp->uncompressed_block_size;
				bytes_written = 0;
				//uint8_t *compressed_buff =  malloc((strlen(char_buff_uncompressed))* sizeof(uint8_t));
				//assert(compressed_buff);	
				//if (rank == master_job_phase_2)
				//	fprintf(stderr, "rank %d :::: start loop compression \n", rank);

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
						// fprintf(stderr, "rank %d :::: block_length = %d \n", rank, block_length);
						//is it necessary?
						//if (block_length < 0) break;

						// count = fwrite(fp->compressed_block, 1, block_length, fp->file);
						// we replace the fwrite with a memcopy
						memcpy(compressed_buff + compressed_size, fp->compressed_block, block_length);
						compressed_size +=(size_t)block_length;
						fp->block_address += block_length;
					}
					//}
				}
			
				free(tmp_buffer);
				//fprintf(stderr, "rank %d :::: finish compression for block = %d compressed size = %zu \n", rank, m6, compressed_size);
				m6++;
			}

				
			//the last rank add the magic number
			//to be compatible with BAM 
			if ( rank == num_proc - 1){
                        	static uint8_t magic[28] =  "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";
				memcpy(compressed_buff + compressed_size, magic, 28);
				compressed_size += 28;
			}	

			BGZF2 *fp_header;
			fp_header = calloc(1, sizeof(BGZF2));
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
				compressed_header =  malloc(strlen(header) * sizeof(uint8_t));

				//fprintf(stderr, "rank %d :::: start loop compression \n", rank);

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
						//fprintf(stderr, "rank %d :::: block_length = %d \n", rank, block_length);
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

			//dispatch header size                               
			MPI_Bcast( &compressed_size_header, 1, MPI_LONG_LONG_INT, master_job_phase_2, COMM_WORLD);
			
			size_t compSize = compressed_size;
		
 			if (rank == master_job_phase_2){
				fprintf(stderr, "Rank %d :::::[WRITE][COMPRESSION] Time for compressing %d blocks in %f seconds \n",
					rank, m6, MPI_Wtime()-time_count);
				fprintf(stderr, "Rank %d :::::[WRITE][COMPRESSION] total compressing size = %zu  \n", rank, compSize);
			}

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
				for (i1 = 1; i1 < (num_proc + 1); i1++) y2[i1] = y[i1-1];
				for (i1 = 1; i1 < (num_proc +1); i1++) y2[i1] = y2[i1-1] + y2[i1];
				for (i1 = 0; i1 < (num_proc +1); i1++) y2[i1] = y2[i1] + write_offset;
				
			}

			MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);

			// BEGIN> FINE TUNING FINFO FOR WRITING OPERATIONS
			//MPI_Info_set(finfo,"striping_factor","12");
			//MPI_Info_set(finfo,"striping_unit","1610612736"); //1G striping
			//MPI_Info_set(finfo,"striping_unit","268435456"); //256 Mo
			//MPI_Info_set(finfo,"nb_proc","12");
			//MPI_Info_set(finfo,"cb_nodes","12");
			//MPI_Info_set(finfo,"cb_block_size","4194304"); /* 4194304 = 4 MBytes - should match FS block size */
			//MPI_Info_set(finfo,"cb_buffer_size","1610612736"); /* 128 MBytes (Optional) */
			// END> FINE TUNING FINFO FOR WRITING OPERATIONS

			ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);

			 MPI_Offset file_size;
                         MPI_File_get_size(out, &file_size);

			if (ierr) {
				fprintf(stderr, "Rank %d failed to open %s.\nAborting.\n\n", rank, path);
				MPI_Abort(COMM_WORLD, ierr);
				exit(2);
			}
			else{
				if(!rank)fprintf(stderr, "Rank %d :::::[WRITE][AFTER COMPRESSION] %s.gz successfully opened\n", rank, chrName);
			}

			time_count = MPI_Wtime();

			if (rank == master_job_phase_2 && (merge == 0)){
				//MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
				MPI_File_write(out, compressed_header, compressed_size_header, MPI_BYTE, MPI_STATUS_IGNORE);
			}
			if (rank == master_job_phase_2 && (merge == 1) && (file_size == 0)){
				// MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
				MPI_File_write(out, compressed_header, compressed_size_header, MPI_BYTE, MPI_STATUS_IGNORE);
			}

			if (!merge) write_offset += compressed_size_header;
                        if  (merge && (file_size == 0)) write_offset += compressed_size_header;
                        else if  (merge && (file_size > 0)) write_offset += file_size;
	
			free(compressed_header);
			time_count=MPI_Wtime();
		
			//we write by block of 1gb
			size_t tmp_size_buffer2 = tmp_size_buffer;

			if (compSize < tmp_size_buffer){

			 	MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
                         	MPI_File_write_all(out, compressed_buff, compSize, MPI_BYTE, &status);
			}
			//we write by block of 1gb
			else {
				uint8_t *compressed_buff_tmp = compressed_buff;
				size_t tmp10 = 0;
				int block_tmp = 0;
				int count_status = 0;
				int error_status = 0;

				while (tmp_size_buffer2 > 0){
				
		 			MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
					MPI_File_write_all(out, compressed_buff_tmp, tmp_size_buffer2, MPI_BYTE, &status);

					MPI_Get_count(&status, MPI_BYTE, &count_status);
					assert(count_status == tmp_size_buffer2);
					compressed_buff_tmp += tmp_size_buffer2;
					tmp10 += tmp_size_buffer2;
					write_offset += tmp_size_buffer2;
					block_tmp++;

					if ( (compSize - tmp10) > tmp_size_buffer ) tmp_size_buffer2 = tmp_size_buffer;
					else tmp_size_buffer2 = (compSize - tmp10);

				}
			}

			MPI_Barrier(COMM_WORLD);
	
			//task FINE TUNING FINFO BACK TO READING OPERATIONS
			//MPI_Info_set(finfo,"striping_factor","12");
			//MPI_Info_set(finfo,"striping_unit","2684354560"); //1G striping

			//MPI_Info_set(finfo,"nb_proc","128");
			//MPI_Info_set(finfo,"cb_nodes","128");
			//MPI_Info_set(finfo,"cb_block_size","2684354560"); /* 4194304 MBytes - should match FS block size */
			//MPI_Info_set(finfo,"cb_buffer_size","2684354560"); /* 128 MBytes (Optional) */
	
			MPI_File_close(&out);	
			
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

			free(compressed_buff);
			free(y);
	                free(y2);
			
		}//end if (write_sam == 1)

		/*
                 * BAM CASE
                 */
		if (write_format == 1){

                        #ifdef HAVE_HTSLIB

				//we make a copy of the header because we modify it
				char *header_tmp = malloc(strlen(header) + 1);
				header_tmp[strlen(header)]=0;
				header_tmp = strcpy(header_tmp, header);

                                int file_exist = 1;
                                if (!merge){
                                        path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
                                        sprintf(path, "%s/%s.bam", output_dir, chrName);
                                        ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);
                                }
                                else{
                                        path = (char*)malloc((strlen(output_dir) + strlen(file_name_sorted) + 40) * sizeof(char));
                                        sprintf(path, "%s/%s", output_dir, file_name_sorted);

                                        file_exist = access(path, F_OK);
                                        ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY|MPI_MODE_CREATE|MPI_MODE_APPEND, finfo, &out);
                                }

                                if (ierr) {
                                        fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] failed to open %s.\nAborting.\n\n", rank, path);
                                        MPI_Abort(COMM_WORLD, ierr);
                                        exit(2);
                                }
                                else{
                                        if((rank == master_job_phase_2) && merge == 0)fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] %s.bam successfully opened\n", rank, chrName);
                                        if((rank == master_job_phase_2) && merge == 1)fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] %s_sorted.bam successfully opened\n", rank, file_name);
                                }


                                sam_hdr_t *h = NULL;
                                hts_idx_t *idx = NULL;
                                bam1_t *b = NULL;
                                int ret = 0;
				int r = 0;	
                                char moder[8];
                                char modew[800];
				uint8_t buff[28];
				hts_opt *in_opts = NULL;
                                hts_opt *out_opts = NULL;
                                struct opts opts;
                                opts.fn_ref = NULL;
                                opts.flag = 0;
                                opts.clevel = -1;
                                opts.ignore_sam_err = 0;
                                opts.nreads = 0;
                                opts.extra_hdr_nuls = 0;
                                opts.benchmark = 0;
				opts.nthreads = 0; 
				opts.multi_reg = 0;
                                opts.index = NULL;
                                opts.min_shift = 0;
                                opts.flag |= WRITE_BINARY_COMP;
                                char buff_tmp[28] = {0};
                                assert(buff_tmp);
                                htsFile *in_header, *in_sam, *out_header_bam, *out_bam;
                                size_t sam_size = strlen(char_buff_uncompressed);
				size_t header_size = strlen(header);
				size_t comp_header_size = 0;
				size_t comp_bam_size = 0;
                                strcpy(moder, "r");
				strcpy(modew, "w");
                                strcat(modew, "b");

				//file descriptor for reading part
                                hFILE *hf_in_header = create_hfile_mem(header_tmp, moder, header_size, header_size);
                                assert(hf_in_header);
                                in_header = hts_hopen(hf_in_header ,"", moder);
				in_header->format.format = sam;	

                                hFILE *hf_in_sam = create_hfile_mem(char_buff_uncompressed, moder, sam_size, sam_size);
                                assert(hf_in_sam);
                                in_sam = hts_hopen(hf_in_sam ,"", moder);
				in_sam->format.format = sam;

				//hfile descriptors for writing part
     				hFILE *hf_out_header = create_hfile_mem(header_tmp, modew, header_size, header_size);
                                assert(hf_out_header);
                                out_header_bam = hts_hopen(hf_out_header ,"", "wb");
				out_header_bam->format.format = bam;


                                hFILE *hf_out_sam = create_hfile_mem(char_buff_uncompressed, modew, sam_size, sam_size);
                                assert(hf_out_sam);
                                out_bam = hts_hopen(hf_out_sam ,"", "wb");
				out_bam->format.format = bam;                                              
 
                                hts_opt_apply(in_header, in_opts);
                                hts_opt_apply(in_sam, in_opts);
				hts_opt_apply(out_header_bam, out_opts);
                                hts_opt_apply(out_bam, out_opts);
	
				//BAM pointer fd
				BGZF *bfp_h = out_header_bam->fp.bgzf;
				BGZF *bfp_b = out_bam->fp.bgzf;

				switch (hts_get_format(in_sam)->category) {
                                	case sequence_data:
                                        	if ( rank == master_job_phase_2)
                                                	fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] we have SAM format in buffer\n", rank);
                                        	break;

                                	default:
                                        	fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] unsupported or unknown category of data in input file\n", rank);

                                }

                                h = sam_hdr_read(in_header);
                                if (h == NULL) {
                                        fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] Couldn't read header for \n", rank); }
                                b = bam_init1();
                                if (b == NULL) {
                                        fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] Out of memory allocating BAM struct \n", rank);}

				MPI_Offset file_size;
                                MPI_File_get_size(out, &file_size);

				//we compress header and write it	
				if ( (rank == master_job_phase_2 && !merge) || (( rank == master_job_phase_2) &&  merge && file_size == 0) ){
                                        if (bam_hdr_write(bfp_h, h) < 0)
                                                fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] Error writing output header.\n", rank);
                                        else
                                                fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] Header commpression ok.\n", rank);

					//add magic in the header bam
					ret = bgzf_raw_write(bfp_h, "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", 28);
					assert(ret == 28);
					//we get the size of compressed header
					char *memptr_h = header_tmp;
                                	while (comp_header_size < ( header_size - 30 )){
                                        	memcpy(buff, memptr_h + comp_header_size , 28);
                                        	if (memcmp(buff,"\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", 28) == 0){
                                                	fprintf(stderr, " Rank %d :::::[WRITE][BAM RESULTS] Header commpression size %zu.\n", rank, comp_header_size);
                                                	break;
                                        	}
                                        	comp_header_size++;
                                	}
                                }

				//dispatch header size                               
				MPI_Bcast( &comp_header_size, 1, MPI_LONG_LONG_INT, master_job_phase_2, COMM_WORLD);                                

                                ret = 0;
                                while ((r = sam_read1(in_sam, h, b)) >= 0)
                                	ret +=bam_write1(bfp_b, b);

				bgzf_flush(bfp_b);
				assert(ret > 0);
         							
                                ret = bgzf_raw_write(bfp_b, "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", 28);
                                assert( ret == 28);
			
				//we compute the length of the compressed SAM buffer
				char *memptr_b = char_buff_uncompressed;
    				while (comp_bam_size < (sam_size - 30)){
        				memcpy(buff, memptr_b + comp_bam_size, 28);
        				if (memcmp(buff,"\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", 28) == 0){
						fprintf(stderr, " Rank %d :::::[WRITE][BAM RESULTS] BAM commpression size %zu.\n", rank, comp_bam_size);                				
                				break;
        				}
       					comp_bam_size++;
    				}
							
				// now the rank can write the buffer hold by memptr
				size_t bamSize = comp_bam_size;

				if ( rank == num_proc - 1)
					bamSize +=28;

                                size_t write_offset = 0;
                                size_t tmp_size_buffer = 1024*1024*1024;
                                size_t tmp_size_buffer2 = tmp_size_buffer;
                                MPI_Offset * y  = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
                                MPI_Offset * y2 = (MPI_Offset *) calloc(num_proc+1, sizeof(MPI_Offset));

                                MPI_Gather(&bamSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);

                                int i1 = 0;
                                if (rank == master_job_phase_2){
                                        for (i1 = 1; i1 < (num_proc + 1); i1++) y2[i1] = y[i1-1];
                                        for (i1 = 1; i1 < (num_proc +1); i1++)  y2[i1] = y2[i1-1] + y2[i1];
                                        for (i1 = 0; i1 < (num_proc +1); i1++)  y2[i1] = y2[i1] + write_offset;

                                }
                                MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);

                                time_count = MPI_Wtime();

				/*
 				 * We write the header 
 				 */  			
				char *memptr_h2 = header_tmp;

				if (rank == master_job_phase_2 && (merge == 0))
					MPI_File_write(out, memptr_h2, comp_header_size, MPI_BYTE, MPI_STATUS_IGNORE);
											
                        	if (rank == master_job_phase_2 && (merge == 1) && (file_size == 0)) 
					MPI_File_write(out, memptr_h2, comp_header_size, MPI_BYTE, MPI_STATUS_IGNORE);
									
				if (!merge) write_offset += comp_header_size;
				if  (merge && (file_size == 0)) write_offset += comp_header_size;
				else if  (merge && (file_size > 0)) write_offset += file_size;

							       
 				//fprintf(stderr, " Rank %d :::::[WRITE][BAM RESULTS] Offset %zu.\n", rank, write_offset);
                                                                                                                        					
				/* 				 
 				 We write the BAM in the output file
 				 */ 
                                MPI_Barrier(COMM_WORLD);
				char *memptr_b2 = char_buff_uncompressed;
                                if (bamSize < tmp_size_buffer){
                                        MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
                                        MPI_File_write_all(out, memptr_b2, bamSize, MPI_BYTE, &status);
                                }

                                //we write by block of 1gb
                                else {
                                        uint8_t *buff_tmp = memptr_b2;
                                        size_t tmp10 = 0;
                                        int block_tmp = 0;
                                        int count_status = 0;
                                        int error_status = 0;
                                        while (tmp_size_buffer2 > 0){
                                                MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
                                                MPI_File_write_all(out, memptr_b2, tmp_size_buffer2, MPI_CHAR, &status);
                                                MPI_Get_count(&status, MPI_CHAR, &count_status);
                                                assert(count_status == tmp_size_buffer2);
                                                buff_tmp += tmp_size_buffer2;
                                                tmp10 += tmp_size_buffer2;
                                                write_offset += tmp_size_buffer2;
                                                block_tmp++;
                                                if ( (bamSize - tmp10) > tmp_size_buffer ) tmp_size_buffer2 = tmp_size_buffer;
                                                else tmp_size_buffer2 = (bamSize - tmp10);
                                        }
                                }


                                if (rank == master_job_phase_2)
                                	MPI_File_close(&out);
				//clean up
				bam_destroy1(b);
				// header_tmp should be freed with sam_hdr_destroy
                                sam_hdr_destroy(h);
				//char_buff_uncompressed shoud be freed with hts_close	
                                ret = hts_close(in_sam);
                                ret = hts_close(in_header);
				/*normally it's enough because other fd point to the same buffer*/

				//free((out_header_bam->fp.bgzf)->uncompressed_block);
				//free(out_bam->fp.bgzf->uncompressed_block);
				bgzf_close(bfp_h);
				bgzf_close(bfp_b);
				free(hf_out_header);
				free(hf_out_sam);
				//ret = hts_close(out_header_bam); //no ok
				//ret = hts_close(out_bam); //no ok
				//ret = hclose(hf_in_sam); //no ok
				//ret = hclose(hf_out_sam); //no ok
				//ret = hclose(hf_out_header); //no ok
                               
                                fprintf(stderr, "Rank %d :::::[WRITE][WRITING BAM] Time for chromosome %s writing %f seconds\n", rank, chrName, MPI_Wtime()-time_count);
 
			#endif
                } //end BAM case 

		//MPI_File_close(&out);
		free(path);
		//if (merge) {free(file_name_tmp);}
		//free(y);
		//free(y2);
		//malloc_trim(0);
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
		int write_format
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

	//Size of the buffer use to compress
	size_t tmp_size_buffer = 1024*1024*1024;
	
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
	//replacement of qksort with stable merge sort
	//qksort(offset_source_index, local_readNum, sizeof(size_t), 0, local_readNum - 1, compare_size_t);
	MergeSortMain( offset_source_index, local_readNum);

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
	
	free(read_size_sorted);
        free(offset_source_index);
        free(offset_source_sorted);

	if ( write_format == 2 ){
		
		time_count = MPI_Wtime();
                size_t write_offset = 0;
                size_t samSize = strlen(char_buff_uncompressed);
                size_t size_header = strlen(header);
                size_t tmp_size_buffer = 1024*1024*1024;
                size_t tmp_size_buffer2 = tmp_size_buffer;

                MPI_Offset * y  = (MPI_Offset *) calloc(split_size, sizeof(MPI_Offset));
                MPI_Offset * y2 = (MPI_Offset *) calloc(split_size+1, sizeof(MPI_Offset));

                MPI_Gather(&samSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, split_comm);

		 int i1 = 0;
                 if ( split_rank == 0 ){
	                 for (i1 = 1; i1 < (split_size + 1); i1++) {
        	                 y2[i1] = y[i1-1];
                         }

                         for (i1 = 1; i1 < (split_size +1); i1++) {
                                 y2[i1] = y2[i1-1] + y2[i1];

                         }
                         for (i1 = 0; i1 < (split_size +1); i1++) {
                                 y2[i1] = y2[i1] + write_offset + size_header;
                         }
                }
                MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, split_comm);

		path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
                sprintf(path, "%s/%s.sam", output_dir, chrName);

		ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);

                if (ierr) {
                fprintf(stderr, "Rank %d failed to open %s.\nAborting.\n\n", split_rank, path);
                	MPI_Abort(COMM_WORLD, ierr);
                        exit(2);
                }
                else{
                	if(!split_rank)fprintf(stderr, "Rank %d :::::[WRITE][OPEN SAM RESULTS] %s.sam successfully opened\n", split_rank, chrName);
                }

                time_count = MPI_Wtime();

                if (split_rank == 0 ) {
                	MPI_File_write(out, header, size_header, MPI_CHAR, MPI_STATUS_IGNORE);
                }
                
		/* we write the SAM */
                if (samSize < tmp_size_buffer){

	                MPI_File_set_view(out, write_offset, MPI_CHAR, MPI_CHAR, "native", finfo);
                        MPI_File_write_all(out, char_buff_uncompressed, samSize, MPI_CHAR, &status);
                }

		//we write by block of 1gb
		else {
			char *buff_tmp = char_buff_uncompressed;
		        size_t tmp10 = 0;
			int block_tmp = 0;
			int count_status = 0;
			int error_status = 0;
	                while (tmp_size_buffer2 > 0){
				MPI_File_set_view(out, write_offset, MPI_CHAR, MPI_CHAR, "native", finfo);
		        	MPI_File_write_all(out, buff_tmp, tmp_size_buffer2, MPI_CHAR, &status);
		        	MPI_Get_count(&status, MPI_CHAR, &count_status);
		        	assert(count_status == tmp_size_buffer2);
		        	buff_tmp += tmp_size_buffer2;
		        	tmp10 += tmp_size_buffer2;
		        	write_offset += tmp_size_buffer2;
		 		block_tmp++;
		        	if ( (samSize - tmp10) > tmp_size_buffer ) tmp_size_buffer2 = tmp_size_buffer;
		        	else tmp_size_buffer2 = (samSize - tmp10);
		 	}
		}

                MPI_File_close(&out);

                if ( split_rank == 0 )
                	fprintf(stderr, "Rank %d :::::[WRITE][WRITING SAM] Time for chromosome %s writing %f seconds\n",
                                        split_rank, chrName, MPI_Wtime()-time_count);

                free(char_buff_uncompressed);
                free(y);
                free(y2);
		free(path);
		
                MPI_Barrier(split_comm);
                
        }
	else {

		char *p6 = char_buff_uncompressed;
		size_t counter_tmp = 0;
		BGZF2 *fp;
        	fp = calloc(1, sizeof(BGZF2));
		uint8_t *compressed_buff =  malloc((strlen(char_buff_uncompressed))* sizeof(uint8_t));
        	assert(compressed_buff);
		int compress_level = compression_level;
		size_t compressed_size = 0;
		time_count = MPI_Wtime();
	
		while (*p6){

            		char *tmp_buffer =malloc(tmp_size_buffer + 1);
			tmp_buffer[tmp_size_buffer]=0;
                	char *p7 = tmp_buffer;
                	size_t counter_tmp = 0;
                	while (*p6 && counter_tmp < tmp_size_buffer) {*p7=*p6; p6++;p7++; counter_tmp++;}

			int block_length = MAX_BLOCK_SIZE;
			int bytes_written;
			//int length = strlen(char_buff_uncompressed);
			int length = counter_tmp;

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

			//const bgzf_byte_t *input = (void *)char_buff_uncompressed;
			//size_t compressed_size = 0;

			if (fp->uncompressed_block == NULL)
	   			fp->uncompressed_block = malloc(fp->uncompressed_block_size);

			//input = (void *)char_buff_uncompressed;
			const bgzf_byte_t *input = (void *)tmp_buffer;
			block_length = fp->uncompressed_block_size;
			bytes_written = 0;
			//char_buff_compressed =  realloc(char_buff_compressed, (strlen(char_buff_uncompressed)+1) * sizeof(uint8_t));
			//assert(char_buff_compressed != 0);
			//char_buff_compressed[strlen(char_buff_uncompressed)]=0;

			while (bytes_written < length) {
				int copy_length = bgzf_min(block_length - fp->block_offset, length - bytes_written);
				bgzf_byte_t* buffer = fp->uncompressed_block;
				memcpy(buffer + fp->block_offset, input, copy_length);
				fp->block_offset += copy_length;
				input += copy_length;
				bytes_written += copy_length;
			
				while (fp->block_offset > 0) {
					int block_length;
					block_length = deflate_block(fp, fp->block_offset);
						//is it necessary?
					//if (block_length < 0) break;
						// count = fwrite(fp->compressed_block, 1, block_length, fp->file);
					// we replace the fwrite with a memcopy
         				memcpy(compressed_buff + compressed_size, fp->compressed_block, block_length);
                                	compressed_size +=(size_t)block_length;
                                	fp->block_address += block_length;
				}			
			} //end data compression
			free(tmp_buffer);
		}
		//we compress the neader
		BGZF2 *fp_header;
		fp_header = calloc(1, sizeof(BGZF2));
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

		//ierr = MPI_File_open(split_comm, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);
		ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);

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
	
		//we write by block of 1gb
		size_t tmp_size_buffer2 = tmp_size_buffer;
		if ( compressed_size < tmp_size_buffer){
	
			MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
			MPI_File_write(out, compressed_buff, (size_t)compressed_size, MPI_BYTE, &status);
		}
		//we write by block of 1gb
		else {
			uint8_t *compressed_buff_tmp = compressed_buff;
	 		size_t tmp10 = 0;
	 		int block_tmp = 0;
	 		int count_status = 0;
	 		int error_status = 0;
	 		while (tmp_size_buffer2 > 0){
	 			time_count=MPI_Wtime();
	        		MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
	        		MPI_File_write_all(out, compressed_buff_tmp, tmp_size_buffer2, MPI_BYTE, &status);
	        		MPI_Get_count(&status, MPI_BYTE, &count_status);
	        		assert(count_status == tmp_size_buffer2);
	        		compressed_buff_tmp += tmp_size_buffer2;
				tmp10 += tmp_size_buffer2;
	        		write_offset += tmp_size_buffer2;
	        		block_tmp++;
	        		if ( (compressed_size - tmp10) > tmp_size_buffer ) tmp_size_buffer2 = tmp_size_buffer;
	        		else tmp_size_buffer2 = compressed_size - tmp10;
	        	
			}
		}

		MPI_Barrier(split_comm);
	
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
		free(compressed_buff);
		free(char_buff_uncompressed);
		free(y);
		free(y2);
		free(path);
		
	} //end if ( write_sam == 1 )

	//malloc_trim(0);
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
		size_t start_offset_in_file,
		int uniq_chr,
		int write_format,
                int merge,
                char file_name_sorted[]
		){


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
	if (uniq_chr) free(data);
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

	/*
	 * for DEBUG
	 *
	fprintf (stderr, "[WRITE] rank = %d :: char_buff_uncompressed = %s \n", rank, char_buff_uncompressed);
	*/
	
	free(new_offset_dest_index_phase3);
	free(data_reads_to_sort);

	
	if (write_format == 2){
		//write format is SAM
		time_count = MPI_Wtime();
		int file_exist = 1;
		size_t write_offset = 0;
                size_t samSize = strlen(char_buff_uncompressed);
                size_t size_header = strlen(header);
                size_t tmp_size_buffer = 1024*1024*1024;
                size_t tmp_size_buffer2 = tmp_size_buffer;
                MPI_Offset * y  = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
                MPI_Offset * y2 = (MPI_Offset *) calloc(num_proc+1, sizeof(MPI_Offset));
                MPI_Gather(&samSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);
		
		int i1 = 0;
                if (rank ==0){
                	for (i1 = 1; i1 < (num_proc + 1); i1++) 
                        	y2[i1] = y[i1-1];
                        

                         for (i1 = 1; i1 < (num_proc +1); i1++) 
                                     y2[i1] = y2[i1-1] + y2[i1];

                         for (i1 = 0; i1 < (num_proc +1); i1++) 
                                     y2[i1] = y2[i1] + write_offset + size_header;
		}                              
                MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);

		if (!merge){
                	path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
                        sprintf(path, "%s/%s.sam", output_dir, chrName);
			ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);
                        }
                else{
	                path = (char*)malloc((strlen(output_dir) + strlen(file_name_sorted) + 40) * sizeof(char));
		        sprintf(path, "%s/%s", output_dir, file_name_sorted);
			file_exist = access(path, F_OK);
                        ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY|MPI_MODE_CREATE|MPI_MODE_APPEND, finfo, &out);
                }
                if (ierr) {
 	               fprintf(stderr, "Rank %d failed to open %s.\nAborting.\n\n", rank, path);
                       MPI_Abort(COMM_WORLD, ierr);
                       exit(2);
                }
                else{
        	        if(!rank && merge == 0)fprintf(stderr, "Rank %d :::::[WRITE][OPEN SAM RESULTS] %s.sam successfully opened\n", rank, chrName);
                        if(!rank && merge == 1)fprintf(stderr, "Rank %d :::::[WRITE][OPEN SAM RESULTS] %s_sorted.sam successfully opened\n", rank, file_name);
                }

                time_count = MPI_Wtime();

                MPI_Offset file_size;
                MPI_File_get_size(out, &file_size);
		if (rank == master_job_phase_2 && (merge == 0)) {
                	MPI_File_write(out, header, size_header, MPI_CHAR, MPI_STATUS_IGNORE);
                }

                if (rank == master_job_phase_2 && (merge == 1) && (file_size == 0)) {
                	MPI_File_write(out, header, size_header, MPI_CHAR, MPI_STATUS_IGNORE);
                }

                //case of merge the file already exist
                if (merge && file_exist == 0) write_offset += (file_size - size_header);
		if (merge && file_exist > 0) write_offset += file_size;

                MPI_Barrier(COMM_WORLD);

		if (samSize < tmp_size_buffer){
                	MPI_File_set_view(out, write_offset, MPI_CHAR, MPI_CHAR, "native", finfo);
                        MPI_File_write_all(out, char_buff_uncompressed, samSize, MPI_CHAR, &status);
                }
                //we write by block of 1gb
                else {
	                char *buff_tmp = char_buff_uncompressed;
                        size_t tmp10 = 0;
                        int block_tmp = 0;
                        int count_status = 0;
                        int error_status = 0;
                        while (tmp_size_buffer2 > 0){
         	               MPI_File_set_view(out, write_offset, MPI_CHAR, MPI_CHAR, "native", finfo);
                               MPI_File_write_all(out, buff_tmp, tmp_size_buffer2, MPI_CHAR, &status);
                               MPI_Get_count(&status, MPI_CHAR, &count_status);
                               assert(count_status == tmp_size_buffer2);
                               buff_tmp += tmp_size_buffer2;
                               tmp10 += tmp_size_buffer2;
                               write_offset += tmp_size_buffer2;
                               block_tmp++;
                               if ( (samSize - tmp10) > tmp_size_buffer ) tmp_size_buffer2 = tmp_size_buffer;
                                        else tmp_size_buffer2 = (samSize - tmp10);
                                }
                 }

                 MPI_File_close(&out);

                 if (rank == master_job_phase_2)
                 	fprintf(stderr, "Rank %d :::::[WRITE][WRITING SAM] Time for chromosome %s writing %f seconds\n",
                        	rank, chrName, MPI_Wtime()-time_count);

                 free(char_buff_uncompressed);
                 free(y);
                 free(y2);
		 free(path);
                 MPI_Barrier(COMM_WORLD);

	}// end write format SAM

	/* 
 	 * BAM case
 	 */	

	if (write_format == 1){

                #ifdef HAVE_HTSLIB
		//we make a copy of the header because we modify it
		char *header_tmp = malloc(strlen(header));
                header_tmp = strcpy(header_tmp, header);

                int file_exist = 1;
                if (!merge){
	        	path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
                        sprintf(path, "%s/%s.bam", output_dir, chrName);
                        ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);
                }
                else{
                        path = (char*)malloc((strlen(output_dir) + strlen(file_name_sorted) + 40) * sizeof(char));
                        sprintf(path, "%s/%s", output_dir, file_name_sorted);

                        file_exist = access(path, F_OK);
                        ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY|MPI_MODE_CREATE|MPI_MODE_APPEND, finfo, &out);
                }
                if (ierr) {
      	        	fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] failed to open %s.\nAborting.\n\n", rank, path);
              	        MPI_Abort(COMM_WORLD, ierr);
                        exit(2);
                }
                else{
	                if((rank == master_job_phase_2) && merge == 0)fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] %s.bam successfully opened\n", rank, chrName);
                        if((rank == master_job_phase_2) && merge == 1)fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] %s_sorted.bam successfully opened\n", rank, file_name);
                }
                sam_hdr_t *h = NULL;
                hts_idx_t *idx = NULL;
                bam1_t *b = NULL;
                int ret = 0;
                int r = 0;
                char moder[8];
                char modew[800];
                uint8_t buff[28];
                hts_opt *in_opts = NULL;
                hts_opt *out_opts = NULL;
                struct opts opts;
                opts.fn_ref = NULL;
                opts.flag = 0;
                opts.clevel = -1;
                opts.ignore_sam_err = 0;
                opts.nreads = 0;
                opts.extra_hdr_nuls = 0;
                opts.benchmark = 0;
                opts.nthreads = 0;
                opts.multi_reg = 0;
                opts.index = NULL;
                opts.min_shift = 0;
                opts.flag |= WRITE_BINARY_COMP;
		char buff_tmp[28] = {0};
                assert(buff_tmp);
                htsFile *in_header, *in_sam, *out_header_bam, *out_bam;
                size_t sam_size = strlen(char_buff_uncompressed);
                size_t header_size = strlen(header);
                size_t comp_header_size = 0;
                size_t comp_bam_size = 0;
                strcpy(moder, "r");
                strcpy(modew, "w");
                strcat(modew, "b");
		//file descriptor for reading part
		hFILE *hf_in_header = create_hfile_mem(header_tmp, moder, header_size, header_size);
                assert(hf_in_header);
                in_header = hts_hopen(hf_in_header ,"", moder);
                in_header->format.format = sam;

                hFILE *hf_in_sam = create_hfile_mem(char_buff_uncompressed, moder, sam_size, sam_size);
                assert(hf_in_sam);
                in_sam = hts_hopen(hf_in_sam ,"", moder);
                in_sam->format.format = sam;

                //hfile descriptors for writing part
               	hFILE *hf_out_header = create_hfile_mem(header_tmp, modew, header_size, header_size);
                assert(hf_out_header);
                out_header_bam = hts_hopen(hf_out_header ,"", "wb");
                out_header_bam->format.format = bam;


                hFILE *hf_out_sam = create_hfile_mem(char_buff_uncompressed, modew, sam_size, sam_size);
                assert(hf_out_sam);
                out_bam = hts_hopen(hf_out_sam ,"", "wb");
                out_bam->format.format = bam;

                hts_opt_apply(in_header, in_opts);
                hts_opt_apply(in_sam, in_opts);
                hts_opt_apply(out_header_bam, out_opts);
                hts_opt_apply(out_bam, out_opts);

                //BAM pointer fd
                BGZF *bfp_h = out_header_bam->fp.bgzf;
                BGZF *bfp_b = out_bam->fp.bgzf;

                switch (hts_get_format(in_sam)->category) {
                        case sequence_data:
                              	if ( rank == master_job_phase_2)
                              	fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] we have SAM format in buffer\n", rank);
                        break;

                        default:
                	fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] unsupported or unknown category of data in input file\n", rank);

                }

                h = sam_hdr_read(in_header);
                if (h == NULL) 
                      	fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] Couldn't read header for \n", rank); 
                b = bam_init1();
                if (b == NULL) 
                        fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] Out of memory allocating BAM struct \n", rank);

                MPI_Offset file_size;
                MPI_File_get_size(out, &file_size);
		//we compress header and write it
		if ( (rank == master_job_phase_2 && !merge) || (( rank == master_job_phase_2) &&  merge && file_size == 0) ){
                        if (bam_hdr_write(bfp_h, h) < 0)
                              	fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] Error writing output header.\n", rank);
                       	else
                       		fprintf(stderr, "Rank %d :::::[WRITE][BAM RESULTS] Header commpression ok.\n", rank);

                        //add magic in the header bam
			ret = bgzf_raw_write(bfp_h, "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", 28);
                        assert(ret == 28);
			//we get the size of compressed header
			char *memptr_h = header_tmp;
                        while (comp_header_size < ( header_size - 30 )){
                              	memcpy(buff, memptr_h + comp_header_size , 28);
                                if (memcmp(buff,"\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", 28) == 0){
                                       	fprintf(stderr, " Rank %d :::::[WRITE][BAM RESULTS] Header commpression size %zu.\n", rank, comp_header_size);
                                        break;
                              	}
                        comp_header_size++;
                	}
                }

		 //dispatch header size                               
		 MPI_Bcast( &comp_header_size, 1, MPI_LONG_LONG_INT, master_job_phase_2, COMM_WORLD);


                ret = 0;
                while ((r = sam_read1(in_sam, h, b)) >= 0)
                    	ret +=bam_write1(bfp_b, b);
                assert(ret > 0);
		bgzf_flush(bfp_b);

                ret = bgzf_raw_write(bfp_b, "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", 28);
                assert( ret == 28);

                //we compute the length of the compressed SAM buffer
                char *memptr_b = char_buff_uncompressed;
                while (comp_bam_size < (sam_size - 30)){
                      	memcpy(buff, memptr_b + comp_bam_size, 28);
                        if (memcmp(buff,"\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", 28) == 0){
                              	fprintf(stderr, " Rank %d :::::[WRITE][BAM RESULTS] BAM commpression size %zu.\n", rank, comp_bam_size);
                                break;
                      	}
                	comp_bam_size++;
                }
                // now the rank can write the buffer hold by memptr
		size_t bamSize = comp_bam_size;
		if ( rank == num_proc - 1)
	                bamSize +=28;

                size_t write_offset = 0;
                size_t tmp_size_buffer = 1024*1024*1024;
                size_t tmp_size_buffer2 = tmp_size_buffer;
                MPI_Offset * y  = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
                MPI_Offset * y2 = (MPI_Offset *) calloc(num_proc+1, sizeof(MPI_Offset));
                MPI_Gather(&bamSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);

                int i1 = 0;
                if (rank == master_job_phase_2){
	        	for (i1 = 1; i1 < (num_proc + 1); i1++) y2[i1] = y[i1-1];
                        for (i1 = 1; i1 < (num_proc +1); i1++)  y2[i1] = y2[i1-1] + y2[i1];
                        for (i1 = 0; i1 < (num_proc +1); i1++)  y2[i1] = y2[i1] + write_offset;

                }
                MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);

                time_count = MPI_Wtime();

                        /*
                         * We write the header 
                    */
		char *memptr_h2 = header_tmp;

                if (rank == master_job_phase_2 && (merge == 0))
                      	MPI_File_write(out, memptr_h2, comp_header_size, MPI_BYTE, MPI_STATUS_IGNORE);

                if (rank == master_job_phase_2 && (merge == 1) && (file_size == 0))
                     	MPI_File_write(out, memptr_h2, comp_header_size, MPI_BYTE, MPI_STATUS_IGNORE);

                if (!merge) write_offset += comp_header_size;
                if  (merge && (file_size == 0)) write_offset += comp_header_size;
                else if  (merge && (file_size > 0)) write_offset += file_size;


		MPI_Barrier(COMM_WORLD);
                //We write the BAM in the output file
		char *memptr_b2 = char_buff_uncompressed;
                if (bamSize < tmp_size_buffer){
                       	MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
                        MPI_File_write_all(out, memptr_b2, bamSize, MPI_BYTE, &status);
                }
		else {
                     	uint8_t *buff_tmp = memptr_b2;
                        size_t tmp10 = 0;
                        int block_tmp = 0;
                        int count_status = 0;
                        int error_status = 0;
                        while (tmp_size_buffer2 > 0){
                             	MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
                                MPI_File_write_all(out, memptr_b2, tmp_size_buffer2, MPI_CHAR, &status);
                                MPI_Get_count(&status, MPI_CHAR, &count_status);
                                assert(count_status == tmp_size_buffer2);
                                buff_tmp += tmp_size_buffer2;
                                tmp10 += tmp_size_buffer2;
                                write_offset += tmp_size_buffer2;
                                block_tmp++;
                                if ( (bamSize - tmp10) > tmp_size_buffer ) tmp_size_buffer2 = tmp_size_buffer;
                	        else tmp_size_buffer2 = (bamSize - tmp10);
                      	}
                }


                if (rank == master_job_phase_2)
                      	MPI_File_close(&out);
		
		bam_destroy1(b);
		sam_hdr_destroy(h);
		ret = hts_close(in_sam);
                ret = hts_close(in_header);

		free(char_buff_uncompressed);
        	        fprintf(stderr, "Rank %d :::::[WRITE][WRITING BAM] Time for chromosome %s writing %f seconds\n", rank, chrName, MPI_Wtime()-time_count);


		free(path);	
	        #endif
        } //end BAM case 

	if (write_format == 0) {

		/*
 		 * BAM case	 
		 */ 

		 int file_exist = 1;
                 if (!merge){
	                path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
                        sprintf(path, "%s/%s.gz", output_dir, chrName);
                        ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);
                 }
                 else{
         		path = (char*)malloc((strlen(output_dir) + strlen(file_name_sorted) + 40) * sizeof(char));
                        sprintf(path, "%s/%s", output_dir, file_name_sorted);

                        file_exist = access(path, F_OK);
                        ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY|MPI_MODE_CREATE|MPI_MODE_APPEND, finfo, &out);
                 }

                 if (ierr) {
                 	fprintf(stderr, "Rank %d :::::[WRITE][BGZF RESULTS] failed to open %s.\nAborting.\n\n", rank, path);
                        MPI_Abort(COMM_WORLD, ierr);
                        exit(2);
                 }
                 else{
                        if((rank == master_job_phase_2) && merge == 0)fprintf(stderr, "Rank %d :::::[WRITE][BGZF RESULTS] %s.gz successfully opened\n", rank, chrName);
                        if((rank == master_job_phase_2) && merge == 1)fprintf(stderr, "Rank %d :::::[WRITE][BGZF RESULTS] %s_sorted.gz successfully opened\n", rank, file_name);
                 }

                 if (rank == master_job_phase_2)
                        fprintf(stderr, "Rank %d :::::[WRITE][BEFORE COMPRESSION] start compression \n", rank);


		time_count = MPI_Wtime();

		BGZF2 *fp;
		fp = calloc(1, sizeof(BGZF2));
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
		assert(compressed_buff);		

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


		//the last rank add the magic number
		//to be compatible with BAM 
		if ( (write_format == 1) && (rank == master_job_phase_2)) {
			static uint8_t magic[28] =  "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";
			memcpy(compressed_buff + compressed_size, magic, 28);
			compressed_size += 28;
		}	

		//we compress the neader

		BGZF2 *fp_header;
		fp_header = calloc(1, sizeof(BGZF2));
		uint8_t *compressed_header = NULL;
		int compressed_size_header = 0;

		if (rank == master_job_phase_2) {

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
			compressed_header =  malloc((strlen(header)+1) * sizeof(uint8_t));
			compressed_header[strlen(header)] = 0;
			
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
		//dispatch header size                               
		MPI_Bcast( &compressed_size_header, 1, MPI_LONG_LONG_INT, master_job_phase_2, COMM_WORLD);
		
		size_t write_offset = 0;

		MPI_Offset * y  = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
		MPI_Offset * y2 = (MPI_Offset *) calloc(num_proc+1, sizeof(MPI_Offset));

		MPI_Gather(&compSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);

		//now we make a cumulative sum
		int i1 = 0;

		if (rank == master_job_phase_2){
			for (i1 = 1; i1 < (num_proc + 1); i1++) y2[i1] = y[i1-1];
			for (i1 = 1; i1 < (num_proc +1); i1++) y2[i1] = y2[i1-1] + y2[i1];
			for (i1 = 0; i1 < (num_proc +1); i1++) y2[i1] = y2[i1] + write_offset;
			
		}

		MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, COMM_WORLD);
		// we create the path where to write for collective write
		

		// BEGIN> FINE TUNING FINFO FOR WRITING OPERATIONS
		/*
		MPI_Info_set(finfo,"striping_factor","128");
		MPI_Info_set(finfo,"striping_unit","1610612736"); //1G striping
		MPI_Info_set(finfo,"nb_proc","128");
		MPI_Info_set(finfo,"cb_nodes","128");
		MPI_Info_set(finfo,"cb_block_size","1610612736"); // 4194304 MBytes - should match FS block size
		MPI_Info_set(finfo,"cb_buffer_size","1610612736"); // 128 MBytes (Optional)
		*/
    		// END> FINE TUNING FINFO FOR WRITING OPERATIONS

		ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);
		MPI_Offset file_size;
                MPI_File_get_size(out, &file_size);
		
		if (ierr) {
			fprintf(stderr, "Rank %d :::::[WRITE_ANY_DIM] failed to open %s.\nAborting.\n\n", rank, path);
			MPI_Abort(COMM_WORLD, ierr);
			exit(2);
		}
		else{
			if(!rank)fprintf(stderr, "Rank %d :::[WRITE_ANY_DIM] %s.bam successfully opened\n", rank, chrName);
		}

		time_count = MPI_Wtime();
		fprintf(stderr, "Rank %d :::::[WRITE_ANY_DIM] we write the header \n", rank);
		if (rank == master_job_phase_2 && (merge == 0)){
			MPI_File_write(out, compressed_header, compressed_size_header, MPI_BYTE, MPI_STATUS_IGNORE);
		}
		if (rank == master_job_phase_2 && (merge == 1) && (file_size == 0)){
			 MPI_File_write(out, compressed_header, compressed_size_header, MPI_BYTE, MPI_STATUS_IGNORE);
		}
		free(compressed_header);		
		
                if (!merge) write_offset += compressed_size_header;
                if  (merge && (file_size == 0)) write_offset += compressed_size_header;
                else if  (merge && (file_size > 0)) write_offset += file_size;
	
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

		if (rank == master_job_phase_2){
			free(fp_header->uncompressed_block);
			free(fp_header->compressed_block);
			free_cache(fp_header);
		}
		free(fp_header);


		MPI_File_close(&out);

		free(path);
		for(m = 0; m < num_proc; m++)
			if (data2[m]) free(data2[m]);
	
		if (data2) free(data2);
		free(offsets_sorted);
		free(new_read_size_sorted_phase3);
		free(data_size_to_sort);
		free(y);
		free(y2);
	}// end write format = bgzf

}





