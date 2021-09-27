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
    sortAnyDim.c

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

#include <sys/mman.h>
#include <sys/stat.h>
#include <assert.h>
#include <ctype.h>
#include <err.h>
#include <fcntl.h>
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <mpi.h>

#include "mergeSort.h"
#include "parser.h"
#include "preWrite.h"
#include "write.h"
#include "mpiSortUtils.h"
#include "writeUtils.h"
#include "qkSort.h"
#include "parallelBitonicSort.h"

void parallel_sort_any_dim(
		int dimensions,
		size_t local_readNum,
		int split_rank,
		int split_size,
		Read **reads,
		int i, //chromosom number
		int chosen_split_rank,
		MPI_Comm split_comm,
		size_t *localReadNumberByChr,
		char *local_data,
		char *file_name,
		char *output_dir,
		MPI_Info finfo,
		int compression_level,
		size_t total_reads_by_chr,
		size_t start_offset_in_file,
		size_t headerSize,
		char* header,
		char *chrNames,
		int mpi_file_split_comm,
		int uniq_chr
		){

	size_t j;

	size_t *local_reads_coordinates_unsorted;
	local_reads_coordinates_unsorted = (size_t*)malloc(local_readNum*sizeof(size_t));
	local_reads_coordinates_unsorted[0] = 0;

	size_t *local_reads_coordinates_sorted;
	local_reads_coordinates_sorted = (size_t*)malloc(local_readNum*sizeof(size_t));
	local_reads_coordinates_sorted[0] = 0;

	int *local_reads_sizes_unsorted; //vector of reads length
	local_reads_sizes_unsorted = (int*)malloc(local_readNum*sizeof(int));
	local_reads_sizes_unsorted[0] = 0;

	int *local_reads_sizes_sorted; //vector of reads length
	local_reads_sizes_sorted = (int*)malloc(local_readNum*sizeof(int));
	local_reads_sizes_sorted[0] = 0;

	int *local_reads_rank_unsorted; //vector of read order in the split_rank
	local_reads_rank_unsorted = (int*)malloc(local_readNum*sizeof(int));
	local_reads_rank_unsorted[0] = 0;

	int *local_reads_rank_sorted; //vector of read order in the split_rank
	local_reads_rank_sorted = (int*)malloc(local_readNum*sizeof(int));
	local_reads_rank_sorted[0] = 0;

	size_t *local_offset_source_unsorted; //vector of read order in the split_rank
	local_offset_source_unsorted = (size_t*)malloc(local_readNum*sizeof(size_t));
	local_offset_source_unsorted[0] = 0;

	size_t *local_offset_source_sorted; //vector of read order in the split_rank
	local_offset_source_sorted = (size_t*)malloc(local_readNum*sizeof(size_t));
	local_offset_source_sorted[0] = 0;

	// the vectors we are going to sort
	// and indexed pbs stands for parallel bitonic sort
	// pbs_vector hold the coordinates and index of the coordinates
	size_t *pbs_local_reads_coordinates;
	size_t *pbs_global_reads_coordinates_index;

	//task Init offset and size for source - free chr
	// from mpiSort_utils.c
	get_coordinates_and_offset_source_and_size_and_free_reads(
			split_rank,
			local_reads_rank_unsorted,
			local_reads_coordinates_unsorted,
			local_offset_source_unsorted,
			local_reads_sizes_unsorted,
			reads[i],
			local_readNum
			);

	//init indices for qksort
	size_t *coord_index = (size_t*)malloc(local_readNum*sizeof(size_t));

	for(j = 0; j < local_readNum; j++){
		coord_index[j] = j;
	}

	//To start we sort locally the reads coordinates.
	//this is to facilitate the bitonic sorting
	//if the local coordinates to sort are to big we could get rid of
	//this step.
	base_arr2 = local_reads_coordinates_unsorted;
	qksort(coord_index, local_readNum, sizeof(size_t), 0, local_readNum - 1, compare_size_t);

	//We index data
	for(j = 0; j < local_readNum; j++){
		local_reads_coordinates_sorted[j] 	= local_reads_coordinates_unsorted[coord_index[j]];
		local_reads_rank_sorted[j] 			= local_reads_rank_unsorted[coord_index[j]];
		local_reads_sizes_sorted[j] 		= local_reads_sizes_unsorted[coord_index[j]];
		local_offset_source_sorted[j] 		= local_offset_source_unsorted[coord_index[j]];
	}

	free(coord_index); //ok
	free(local_reads_rank_unsorted); //ok
	free(local_reads_coordinates_unsorted); //ok
	free(local_reads_sizes_unsorted); //ok
	free(local_offset_source_unsorted); //ok

	MPI_Barrier(split_comm); //maybe not necessary

	//first we get the vector of all offset in destination file
	size_t *all_reads_coordinates 			= NULL;
	size_t *all_reads_coordinates_sorted 	= NULL;
	size_t *all_offsets_sources 			= NULL;
	int *all_reads_sizes 					= NULL;
	int *all_reads_rank 					= NULL;
	size_t *all_reads_coordinates_index 	= NULL;

	//we get the total number of reads.
	//because the elected is going to collect them
	//and dispatch to ranks belonging to the bitonic dimension
	size_t total_num_read = 0;
	MPI_Reduce(
			&localReadNumberByChr[i],
			&total_num_read,
			1,
			MPI_LONG_LONG_INT,
			MPI_SUM,
			chosen_split_rank,
			split_comm
			);

	if (split_rank == chosen_split_rank)
			fprintf(stderr,	"rank %d :::::[SORT_ANY_DIM] total_num_read = %zu \n", split_rank, total_num_read);


	MPI_Barrier(split_comm); //maybe not necessary

	if (split_rank == chosen_split_rank){
		all_reads_coordinates 	= malloc (total_num_read * sizeof(size_t));
		all_offsets_sources 	= malloc (total_num_read * sizeof(size_t));
		all_reads_sizes 		= malloc (total_num_read * sizeof(int));
		all_reads_rank 			= malloc (total_num_read * sizeof(int));

		size_t k =0;

		for (k = 0; k < total_num_read; k++){
			all_reads_sizes[k] 			= 0;
			all_reads_coordinates[k] 	= 0;
			all_offsets_sources[k] 		= 0;
			all_reads_rank[k] 			= 0;
		}
	}


	/*
	 *
	 *
	 * But each jobs has vector of diffrent length to sort.
	 *
	 */


	MPI_Barrier(split_comm); //maybe not necessary

	// we broadcast the total number of reads to each rank
	MPI_Bcast(&total_num_read, 1, MPI_LONG_LONG_INT, chosen_split_rank, split_comm );


	// vector of number of read per jobs
	size_t *num_reads_per_jobs = (size_t *) malloc(split_size * sizeof(size_t));
	// start_num_reads_per_jobs is the start index when dispatching the total
	size_t *start_num_reads_per_jobs = (size_t *) malloc((split_size + 1) * sizeof(size_t));

	// Preparation of recieving the coordinates, offset, and reads size.
	// chosen_rank recieves the number of reads of each rank and put it in a vector
	MPI_Gather(
			&local_readNum,
			1,
			MPI_LONG_LONG_INT,
			&num_reads_per_jobs[split_rank - chosen_split_rank],
			1,
			MPI_LONG_LONG_INT,
			chosen_split_rank,
			split_comm
			);

	MPI_Barrier(split_comm);

	// we initialyze start_num_reads_per_jobs
	if (split_rank == chosen_split_rank){

		start_num_reads_per_jobs[0] = 0;
		int k = 0;
		for (k = 1; k < (split_size + 1); k++){
			start_num_reads_per_jobs[k] = num_reads_per_jobs[k-1];
		}

		for (k = 1; k < split_size; k++){
			size_t tmp = start_num_reads_per_jobs[k - 1];
			size_t tmp2 = start_num_reads_per_jobs[k];
			start_num_reads_per_jobs[k] = tmp + tmp2;
		}
	}

	/*
	 * split_chosen_rank collects sizes, coordinates and offsets
	 * in all_vector_
	 */
	double time_count  = MPI_Wtime();
	double time_count1 = MPI_Wtime();

	if (split_rank ==chosen_split_rank){

		MPI_Status status;
		//we copy the first elements in
		int k=0;
		size_t st = start_num_reads_per_jobs[chosen_split_rank];
		for (k = 0; k < num_reads_per_jobs[chosen_split_rank]; k++){
			all_reads_rank[st] 			= local_reads_rank_sorted[k];
			all_reads_sizes[st] 		= local_reads_sizes_sorted[k];
			all_offsets_sources[st] 	= local_offset_source_sorted[k];
			all_reads_coordinates[st] 	= local_reads_coordinates_sorted[k];
			st++;
		}

		for(j = 0; j < split_size; j++){
			if(j != chosen_split_rank){

				// first we care for ranks
				int *temp_buf 		= malloc(num_reads_per_jobs[j]* sizeof(int));
				int *temp_buf1 		= malloc(num_reads_per_jobs[j]* sizeof(int));
				size_t *temp_buf2 	= malloc(num_reads_per_jobs[j]* sizeof(size_t));
				size_t *temp_buf3 	= malloc(num_reads_per_jobs[j]* sizeof(size_t));

				MPI_Recv(temp_buf, num_reads_per_jobs[j], MPI_INT, j, 0, split_comm, &status);
				MPI_Recv(temp_buf1, num_reads_per_jobs[j], MPI_INT, j, 1, split_comm, &status);
				MPI_Recv(temp_buf2, num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 2, split_comm, &status);
				MPI_Recv(temp_buf3, num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 3, split_comm, &status);

				st=0;
				size_t st = start_num_reads_per_jobs[j];

				for (k = 0; k < num_reads_per_jobs[j]; k++){
					all_reads_rank[st] = temp_buf[k];
					all_reads_sizes[st] = temp_buf1[k];
					all_offsets_sources[st] = temp_buf2[k];
					all_reads_coordinates[st] = temp_buf3[k];
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
		MPI_Send(local_reads_rank_sorted, local_readNum, MPI_INT, chosen_split_rank,  0, split_comm);
		MPI_Send(local_reads_sizes_sorted, local_readNum, MPI_INT, chosen_split_rank,  1, split_comm);
		MPI_Send(local_offset_source_sorted, local_readNum, MPI_LONG_LONG_INT, chosen_split_rank,  2, split_comm);
		MPI_Send(local_reads_coordinates_sorted, local_readNum, MPI_LONG_LONG_INT, chosen_split_rank,  3, split_comm);

	}

	if (split_rank == chosen_split_rank)
		fprintf(stderr,	"rank %d :::::[SORT_ANY_DIM] Time to collect before bitonic part  %f seconds\n", split_rank, MPI_Wtime() - time_count);

	free(local_reads_sizes_sorted); //ok
	free(local_reads_rank_sorted); //ok
	free(local_offset_source_sorted); //ok
	free(local_reads_coordinates_sorted); //ok

	/*
	 * ENTER BITONIC PART
	 *
	 * The reads are sorted according to their coordinates.
	 *
	 * Now rank 0 hold all the coordinates, offset and size.
	 *
	 * The bitonic works with a power of 2 number ranks and works with vectors of equal lentgh.
	 *
	 * First we compute a dimensions (2^n) and a length of vectors before dispaching them
	 * to the bitonic rank.
	 *
	 * Improvement:
	 *
	 * This part is not optimal because a rank is going to recieve
	 * all the offset.
	 *
	 * What we should do is to equilibrate the vector by trading
	 * data between ranks in a ring fashion. This should eliminate the pressure
	 * on elected rank.
	 *
	 */

	// we compute the length of the vector to recieve
	// the we split all_offset_dest_index_phase1 among the dimension processors
	size_t pbs_num_coordinates_to_recieve = total_num_read/dimensions;
	// for the last processor will recieve the rest of the division
	size_t pbs_num_coordinates_to_recieve_left  =  total_num_read - pbs_num_coordinates_to_recieve * dimensions;


	/*
	 * Here we reduce the dimension.
	 * We do that because we don't want a job to finish with zero vector
	 * after the bitonic.
	 *
	 * In the future we shall use a bruck to dispatch the read between dimensions
	 */

	while((pbs_num_coordinates_to_recieve_left * dimensions) > pbs_num_coordinates_to_recieve){
		dimensions >>= 1;
		pbs_num_coordinates_to_recieve = total_num_read/dimensions;
		pbs_num_coordinates_to_recieve_left = total_num_read - pbs_num_coordinates_to_recieve * dimensions;

	}

	// we compute a vector of size dimensions which contain the number
	// of reads to send
	size_t *pbs_local_num_read_per_job = (size_t *)malloc(dimensions * sizeof(size_t));

	// each job create a vector with with the length
	// of the offset vector for each job in [0, dimension[
	for (j = 0; j < dimensions ; j++){
		// we add pbs_num_coordinates_to_recieve_left
		// because we want the vector to have the size size
		pbs_local_num_read_per_job[j] = pbs_num_coordinates_to_recieve + pbs_num_coordinates_to_recieve_left;
	}

	// the lastest rank get the reads that left
	size_t *pbs_start_num_coordinates_per_jobs = (size_t *) malloc((dimensions + 1) * sizeof(size_t));

	// the master job compute the start index of the element
	// to dispatch
	if (split_rank == chosen_split_rank){

		pbs_start_num_coordinates_per_jobs[0] = 0;
		int k = 0;
		for (k = 1; k < (dimensions +1); k++){
			pbs_start_num_coordinates_per_jobs[k] = pbs_local_num_read_per_job[k-1];
		}
		for (k = 1; k < dimensions; k++){
			size_t tmp = pbs_start_num_coordinates_per_jobs[k - 1];
			size_t tmp2 = pbs_start_num_coordinates_per_jobs[k];
			// we remove the left over reads
			pbs_start_num_coordinates_per_jobs[k] = tmp + tmp2 - pbs_num_coordinates_to_recieve_left;
		}
	}

	// the processors chosen rank send the coordinates
	// and read sizes to all the rank in [0-dimension]

	MPI_Barrier(split_comm);

	if (split_rank < dimensions){

		// pbs_local_reads_coordinates is a table containing the unsorted
		// reference coordinates
		pbs_local_reads_coordinates = (size_t *)malloc(sizeof(size_t) * pbs_local_num_read_per_job[split_rank]);
		pbs_local_reads_coordinates[0] = 0;

		//now the master send
		time_count = MPI_Wtime();
		if ( split_rank != chosen_split_rank ){

				MPI_Status status;
				MPI_Recv(
						pbs_local_reads_coordinates,
						pbs_local_num_read_per_job[split_rank],
						MPI_LONG_LONG_INT,
						chosen_split_rank,
						0,
						split_comm,
						&status
						);

		}
		else {
			//first we copy the data from the master job
			size_t ind = pbs_start_num_coordinates_per_jobs[chosen_split_rank];
			int k = 0;

			for (k = 0; k < pbs_local_num_read_per_job[chosen_split_rank]; k++){
				pbs_local_reads_coordinates[k] = all_reads_coordinates[ind];

				ind++;
			}

			for(j = 0; j < dimensions; j++){

				if (j != chosen_split_rank){

					MPI_Send(
							&all_reads_coordinates[pbs_start_num_coordinates_per_jobs[j]],
							pbs_local_num_read_per_job[j],
							MPI_LONG_LONG_INT,
							j,
							chosen_split_rank,
							split_comm
							);
				}
			}
		}

		if (split_rank == chosen_split_rank){
			free(all_reads_coordinates);
		}


		// we build pbs_local_reads_coordinates_index
		// the index we be used to reoder offset and read size
		// to comput the offset destination
		pbs_global_reads_coordinates_index = (size_t *)malloc(pbs_local_num_read_per_job[split_rank]*sizeof(size_t));
		pbs_global_reads_coordinates_index[0] = 0;

		for (j = 0; j < pbs_local_num_read_per_job[split_rank]; j++){

			if (split_rank == chosen_split_rank){
				pbs_global_reads_coordinates_index[j] = j + pbs_local_num_read_per_job[split_rank] * split_rank;
			}
			else{
				pbs_global_reads_coordinates_index[j] = j + pbs_local_num_read_per_job[split_rank] * split_rank -
						(split_rank * pbs_num_coordinates_to_recieve_left);
			}
		}

		// now each rank from [0, dimension[
		// is going to bitonic sort
		// input are:
		// pbs_local_reads_coordinates
		// pbs_local_reads_coordinates_index
		time_count = MPI_Wtime();

		if (split_rank == chosen_split_rank)
			fprintf(stderr,	"rank %d :::::[SORT_ANY_DIM] Call bitonic for reads coordinates with dimensions = %d \n", split_rank, dimensions);

		ParallelBitonicSort(
				split_comm,
				split_rank,
				dimensions,
				pbs_local_reads_coordinates,
				pbs_global_reads_coordinates_index,
				pbs_local_num_read_per_job[split_rank],
				pbs_num_coordinates_to_recieve_left
				);

		size_t k1 = 0;
		for (k1 = 1; k1 < pbs_local_num_read_per_job[split_rank]; k1++){
			assert( pbs_local_reads_coordinates[k1 - 1] <= pbs_local_reads_coordinates[k1]);
		}

		if (split_rank == chosen_split_rank)
			fprintf(stderr,	"rank %d :::::[SORT_ANY_DIM] Time in parallel bitonic sort  %f seconds\n", split_rank, MPI_Wtime() - time_count);

		//we compute a new total number of reads
		size_t total_num_read_after_bitonic_sort = 0;
		int k=0;
		for (k = 0; k < dimensions; k++){
			total_num_read_after_bitonic_sort += pbs_local_num_read_per_job[k];
		}

		// now we gather all the pbs_local_dest_offset_index
		// and pbs_local_dest_offset in 2 vectors
		// all_offset_dest_sorted_phase1
		// all_offset_index_phase_1

		//we allocate vector to send
		// we remove zero
		size_t start_index=0;

		while (pbs_local_reads_coordinates[start_index] == 0){
			start_index++;
		}

		pbs_local_num_read_per_job[split_rank] -=  start_index;
		all_reads_coordinates_index 			= (size_t *)malloc(sizeof(size_t) * total_num_read);
		all_reads_coordinates_sorted 			= (size_t *)malloc (total_num_read * sizeof(size_t));

		if (split_rank == chosen_split_rank){

			pbs_start_num_coordinates_per_jobs[0] = 0;

			for (k = 1; k < (dimensions +1); k++){
				pbs_start_num_coordinates_per_jobs[k] = pbs_local_num_read_per_job[k-1];
			}
			for (k = 1; k < dimensions; k++){
				size_t tmp = pbs_start_num_coordinates_per_jobs[k - 1];
				size_t tmp2 = pbs_start_num_coordinates_per_jobs[k];
				// we remove the left over reads
				pbs_start_num_coordinates_per_jobs[k] = tmp + tmp2;
			}
		}


		time_count = MPI_Wtime();

		chosen_split_rank_gather_size_t(
				split_comm,
				split_rank, dimensions, chosen_split_rank,
				pbs_local_num_read_per_job[split_rank],
				pbs_local_num_read_per_job,
				pbs_start_num_coordinates_per_jobs,
				all_reads_coordinates_index,
				pbs_global_reads_coordinates_index,
				start_index
				);

		chosen_split_rank_gather_size_t(
				split_comm,
				split_rank,
				dimensions,
				chosen_split_rank,
				pbs_local_num_read_per_job[split_rank],
				pbs_local_num_read_per_job,
				pbs_start_num_coordinates_per_jobs,
				all_reads_coordinates_sorted,
				pbs_local_reads_coordinates,
				start_index
				);


		free(pbs_local_reads_coordinates); //ok


		//now we verify that all_reads_coordinates_sorted is sorted
		if (split_rank == chosen_split_rank){
			for (k = 1; k < total_num_read ; k++){
				assert(all_reads_coordinates_sorted[k-1] <= all_reads_coordinates_sorted[k]);
			}
		}
		if (split_rank == chosen_split_rank)
			fprintf(stderr,	"rank %d :::::[SORT_ANY_DIM] Time to gather all the reads coordinates %f seconds\n",
					split_rank, MPI_Wtime() - time_count);

		free(pbs_global_reads_coordinates_index); //ok


	} //end if (split_rank < dimensions)

	free(pbs_local_num_read_per_job);
	free(pbs_start_num_coordinates_per_jobs); //ok

	/*
	 * Bitonic is finished
	 *
	 * Now we re-order and dispatch
	 * the ranks, source offset and sizes
	 *
	 * Improvement.
	 *
	 * This step is not necessary. The rordering of the offset sources, rank
	 * and read sizes could be done during the bitonic.
	 *
	 */


	// we create a datatype by split_rank
	size_t *all_offset_dest_sorted		= NULL;
	int *all_reads_size_sorted 			= NULL;
	int *all_reads_rank_sorted 			= NULL;
	size_t *all_offsets_sources_sorted 	= NULL;

	if (split_rank == chosen_split_rank){

		size_t k;
		//we reorder all_reads_sizes
		all_reads_size_sorted 		= malloc(sizeof(int) * total_num_read);
		all_reads_rank_sorted		= malloc(sizeof(int) * total_num_read);
		all_offsets_sources_sorted 	= malloc(sizeof(size_t) * total_num_read);

		for (k = 0; k < total_num_read ; k++){
			all_reads_size_sorted[k] 		= all_reads_sizes[all_reads_coordinates_index[k]];
			all_offsets_sources_sorted[k] 	= all_offsets_sources[all_reads_coordinates_index[k]];
			all_reads_rank_sorted[k] 		= all_reads_rank[all_reads_coordinates_index[k]];

		}

		free(all_reads_sizes); //ok
		free(all_offsets_sources); //ok
		free(all_reads_rank); //ok

		// all_offset_dest holds the output offset of sorted reads
		all_offset_dest_sorted 		= malloc(sizeof(size_t) * total_num_read);
		all_offset_dest_sorted[0] 	= all_reads_size_sorted[0];
		all_offset_dest_sorted[0]  += headerSize;

		for (k = 1; k < total_num_read ; k++){
			all_offset_dest_sorted[k] = all_offset_dest_sorted[k - 1] + all_reads_size_sorted[k];
		}


	} //end if (split_rank == chosen_rank)

	size_t *local_dest_offsets_sorted 		= malloc(sizeof(size_t)*local_readNum);
	size_t *local_source_offsets_sorted 	= malloc(sizeof(size_t)*local_readNum);
	int *local_read_size_sorted 			= malloc(sizeof(size_t)*local_readNum);
	int *local_rank_sorted 					= malloc(sizeof(size_t)*local_readNum);

	/*
	 * now we dispatch all_offset_dest
	 * all_offset_source
	 * all_reads_size
	 */

	time_count = MPI_Wtime();

	if (split_rank != chosen_split_rank){
		//fprintf(stderr, "%d ::::: [send_size_t_master_to_all] rank %d recv %zu from %d \n",rank, rank, size, master);
		MPI_Recv(local_dest_offsets_sorted, local_readNum, MPI_LONG_LONG_INT, chosen_split_rank, 0, split_comm, MPI_STATUS_IGNORE);
		MPI_Recv(local_source_offsets_sorted, local_readNum, MPI_LONG_LONG_INT, chosen_split_rank, 1, split_comm, MPI_STATUS_IGNORE);
		MPI_Recv(local_read_size_sorted, local_readNum, MPI_INT, chosen_split_rank, 2, split_comm, MPI_STATUS_IGNORE);
		MPI_Recv(local_rank_sorted, local_readNum, MPI_INT, chosen_split_rank, 3, split_comm, MPI_STATUS_IGNORE);
	}
	else {
		size_t k=0;
		size_t ind = start_num_reads_per_jobs[chosen_split_rank];

		for (k = 0; k < (num_reads_per_jobs[chosen_split_rank]); k++){

			local_dest_offsets_sorted[k] 	= all_offset_dest_sorted[ind];
			local_source_offsets_sorted[k] 	= all_offsets_sources_sorted[ind];
			local_read_size_sorted[k] 		= all_reads_size_sorted[ind];
			local_rank_sorted[k] 			= all_reads_rank_sorted[ind];

			ind++;
		}

		for(j = 0; j < split_size; j++){

			if (j != chosen_split_rank){
				//fprintf(stderr, "%d ::::: [send_size_t_master_to_all] rank %d send %zu to %d from %zu\n",
				//	rank, rank, size_per_jobs[j], j, start_size_per_job[j]);
				MPI_Send(&all_offset_dest_sorted[start_num_reads_per_jobs[j]],
						num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, split_comm);

				MPI_Send(&all_offsets_sources_sorted[start_num_reads_per_jobs[j]],
						num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 1, split_comm);

				MPI_Send(&all_reads_size_sorted[start_num_reads_per_jobs[j]],
						num_reads_per_jobs[j], MPI_INT, j, 2, split_comm);

				MPI_Send(&all_reads_rank_sorted[start_num_reads_per_jobs[j]],
						num_reads_per_jobs[j], MPI_INT, j, 3, split_comm);

			}
		}
	}

	if (split_rank == chosen_split_rank)
		fprintf(stderr,	"rank %d :::::[SORT_ANY_DIM] Time dispatch all destination offsets %f seconds\n", split_rank, MPI_Wtime() - time_count);

	free(start_num_reads_per_jobs); // ok
	free(num_reads_per_jobs); // ok


	if (split_rank == chosen_split_rank){
		free(all_reads_rank_sorted);
		free(all_offsets_sources_sorted);
		free(all_offset_dest_sorted);
		free(all_reads_size_sorted);
	}


	MPI_Barrier(split_comm);

	if (split_rank < dimensions){
		free(all_reads_coordinates_index); //ok
		free(all_reads_coordinates_sorted);
	}

	if (split_rank == chosen_split_rank)
		fprintf(stderr,	"rank %d :::::[SORT_ANY_DIM] CALL WRITE SAM ANY DIMENSION\n", split_rank);

	writeSam_any_dim(
			dimensions,
			split_rank,
			output_dir,
			header,
			local_readNum,
			total_reads_by_chr,
			chrNames,
			reads[i],
			split_size,
			split_comm,
			chosen_split_rank,
			file_name,
			mpi_file_split_comm,
			finfo,
			compression_level,
			local_dest_offsets_sorted,
			local_source_offsets_sorted,
			local_read_size_sorted,
			local_rank_sorted,
			local_data,
			start_offset_in_file,
			uniq_chr
			);

}
