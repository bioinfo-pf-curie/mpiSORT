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
     write2.c

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/




#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "mpi.h"
#include "parser.h"
#include "assert.h"
#include "preWrite.h"
#include "write2.h"
#include "diffuse.h"
#include "math.h"
#include "sys/mman.h"
#include "bgzf.c"
#include "bgzf.h"
//#include "parabitonicsort.h"

#define MBYTE 1048576

size_t *base_arr2;

static int compare_size_t(const void *a, const void *b){

	 size_t aa = *(size_t *)a, bb = *(size_t *)b;

	 if (base_arr2[aa] > base_arr2[bb])
		return 1;
	else if (base_arr2[aa] < base_arr2[bb])
		return -1;
	else
		return 0;

}

static int compare_size_t_V2(const void *a, const void *b){

	if (*(const size_t *)a > *(const size_t *)b)
		return 1;
	else if (*(const size_t *)a < *(const size_t *)b)
		return -1;
	else
		return 0;

}

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
		int *new_size, int ***data_size
	)
{

	bruck_reads(rank, num_proc, buffs_by_procs, *data2);
	bruck_size(rank, num_proc, local_readNum, number_of_reads_by_procs, *data_size, new_rank, new_size);
	bruck_offsets(rank, num_proc, local_readNum, number_of_reads_by_procs, *data_offsets, new_rank, new_offset);

}

void bruck_reads(int rank, int num_proc, size_t * buffs_by_procs, char** data2)
{
	MPI_Comm comm = MPI_COMM_WORLD;

	int k, m, srank, rrank;
	MPI_Datatype dt_send, dt_recv;
	size_t *recv_size_by_proc=NULL, *send_size_by_proc=NULL;
	int *recv_index=NULL;
	size_t total, send_total;
	int packsize;
	double time;

	int count;

	time = MPI_Wtime();
	for(k=1; k<num_proc; k<<=1)
	{
		srank = (rank - k + num_proc) % num_proc;	//Rank to send to
		rrank = (rank + k) % num_proc;	//Rank to recv from

		count = create_send_datatype_for_reads(rank, num_proc, buffs_by_procs, data2, k, &dt_send, &recv_index);
		MPI_Type_commit(&dt_send);

		send_size_by_proc = (size_t*)malloc(count*sizeof(size_t));
		recv_size_by_proc = (size_t*)malloc(count*sizeof(size_t));

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

		char* interbuff = (char*)calloc(total, sizeof(char));

		MPI_Sendrecv(MPI_BOTTOM, 1, dt_send, srank, 0,
				interbuff, total, MPI_PACKED, rrank, 0, comm, MPI_STATUS_IGNORE);

		for ( m = 0; m < count; m++){
			// we free and allocate data2
			// according to the recieve size
			if (data2[recv_index[m]]){
				free(data2[recv_index[m]]);
				data2[recv_index[m]] = (char *)malloc(sizeof(char)*(recv_size_by_proc[m])+1);
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
			buffs_by_procs[recv_index[m]] = strlen(data2[recv_index[m]]);
		}

		MPI_Barrier(comm);
		MPI_Type_free(&dt_recv);
		MPI_Type_free(&dt_send);

		count = 0;
		free(interbuff);
		free(recv_index);
		free(recv_size_by_proc);
		free(send_size_by_proc);
	}

	/*
	fprintf(stderr, "rank %d ::::: [BRUCK READS] data2 = \n", rank);
		for (m = 0; m < buffs_by_procs[rank]; m++){
			fprintf(stderr, "%c", data2[rank][m]);
		}
	*/
}

void bruck_offsets(int rank, int num_proc, int local_readNum, size_t* number_of_reads_by_procs, size_t ** data_offsets, int *new_rank, size_t* new_offset)
{
	MPI_Comm comm = MPI_COMM_WORLD;

	int k, m, j, srank, rrank;
	MPI_Datatype dt_send;

	size_t *recv_size_by_proc=NULL, *send_size_by_proc=NULL;
	int *recv_index=NULL;
	size_t total, send_total;
	int packsize;
	double time;
	MPI_Status status;

	int count;

	time = MPI_Wtime();

	for (m = 0; m < num_proc; m++){
		number_of_reads_by_procs[m] = 0;

	}

	for(m = 0; m < local_readNum; m++){
		number_of_reads_by_procs[new_rank[m]]++;
	}

	size_t **data_offsets2 = (size_t **)malloc(sizeof(size_t *)*num_proc);
	//we initialize data_offsets
	for(m = 0; m < num_proc; m++){
		data_offsets[m] = NULL;
		data_offsets2[m] = calloc( number_of_reads_by_procs[m], sizeof(size_t));
	}

	size_t *read_by_proc_offset = (size_t *)calloc(num_proc, sizeof(size_t));

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

		//fprintf(stderr, "Rank %d ::::: recv from = %d ::: send to  = %d \n", rank, srank, rrank);

		count = create_send_datatype_for_offsets(rank, num_proc, number_of_reads_by_procs,
				data_offsets, k, &dt_send, &recv_index);

		MPI_Type_commit(&dt_send);

		send_size_by_proc = (size_t*)malloc(count*sizeof(size_t));
		recv_size_by_proc = (size_t*)malloc(count*sizeof(size_t));

		send_total = get_send_size(rank, num_proc, number_of_reads_by_procs,
				&send_size_by_proc, count, k);

		MPI_Pack_size(1, dt_send, comm, &packsize);
		/*
		fprintf(stderr, "Rank %d ::::: IN BRUCK OFFSET ::: packsize = %d ::: send_total = %zu :::: count = %d \n", rank,
				packsize, send_total, count);
		*/

		assert(packsize == (8 * send_total));

		MPI_Sendrecv(send_size_by_proc, count, MPI_LONG_LONG_INT, srank, 0,
			recv_size_by_proc, count, MPI_LONG_LONG_INT,
				rrank, 0, comm, MPI_STATUS_IGNORE);

		total = 0;

		for(m = 0; m < count; m++)
		{
			total += recv_size_by_proc[m];
		}

		//fprintf(stderr, "Rank %d ::::: IN BRUCK OFFSET AFTER SEND RECV ::: total = %zu\n", rank, total);

		size_t *interbuff_offset = calloc(total, sizeof(size_t));
		//interbuff_offset = realloc(interbuff_offset, total);
		MPI_Sendrecv(MPI_BOTTOM, 1, dt_send, srank, 0,
				interbuff_offset, total, MPI_LONG_LONG_INT, rrank, 0, comm, &status);
		/*
		for ( m = 0; m < total; m++){
			fprintf(stderr, "Rank %d ::::: interbuff_offset[%d] = %zu\n", rank, m, interbuff_offset[m]);
		}

		int count_ptr = 0;
		MPI_Get_count(&status, MPI_INT, &count_ptr);

		fprintf(stderr, "Rank %d ::::: IN BRUCK OFFSET ::: statuc count = %d\n", rank, count_ptr);

		if (rank == 2){
			for ( m = 0; m < count; m++){
				fprintf(stderr, "Rank %d ::::: IN BRUCK OFFSET ::: recv_index[%d] = %d\n", rank, m, recv_index[m]);
			}
		}
		 */
		for ( m = 0; m < count; m++){
			// we free and allocate data_offsets
			// according to the recieve size
			if (data_offsets[recv_index[m]]){

				free(data_offsets[recv_index[m]]);

				data_offsets[recv_index[m]] = NULL;

				if (recv_size_by_proc[m] > 0)
					data_offsets[recv_index[m]] = (size_t *)malloc(sizeof(size_t)*(recv_size_by_proc[m]));

				//data_offsets[recv_index[m]] = realloc(data_offsets[recv_index[m]], sizeof(size_t)*(recv_size_by_proc[m]));
				if (data_offsets[recv_index[m]])
					data_offsets[recv_index[m]][0] = 0;
			}
		}
		/*
		MPI_Aint indices_offset[count];
		int blocklens_offset[count];
		MPI_Datatype oldtypes_offset[count];

		for (m = 0; m < count; m++){
			blocklens_offset[m] = (int)(recv_size_by_proc[m]);
			MPI_Get_address(data_offsets[recv_index[m]], &indices_offset[m]);
			oldtypes_offset[m] = MPI_LONG_LONG_INT;
			number_of_reads_by_procs[recv_index[m]] = recv_size_by_proc[m];
		}

		//Create structure of recieve type
		MPI_Type_create_struct(1, blocklens_offset, indices_offset, oldtypes_offset, &dt_recv);
		MPI_Type_commit(&dt_recv);

		int pos=0;
		MPI_Unpack(interbuff_offset, total, &pos, MPI_BOTTOM, 1, dt_recv, comm);
		*/

		size_t *tmp_var = interbuff_offset;

		for (m = 0; m < count; m++){

			memcpy(data_offsets[recv_index[m]], tmp_var, recv_size_by_proc[m] * sizeof(size_t));
			tmp_var += recv_size_by_proc[m];
			number_of_reads_by_procs[recv_index[m]] = recv_size_by_proc[m];

		}
		/*
		fprintf(stderr, "Rank %d ::::: Rank 2 recieve from %d \n", rank, rrank);
		fprintf(stderr, "Rank %d ::::: Rank 2 send to %d \n", rank, srank);
		for (m = 0; m < count; m++){
			for (j = 0; j< number_of_reads_by_procs[recv_index[m]]; j++){
				if (data_offsets[recv_index[m]][j] == 0){
					fprintf(stderr, "Rank %d ::::: data_offsets[recv_index[%d]][%d] = 0 \n", rank, m, j);
				}
			}
		}
		*/

		for (m = 0; m < count; m++){
			for (j = 0; j< number_of_reads_by_procs[recv_index[m]]; j++){
				assert(data_offsets[recv_index[m]][j] != 0);
			}
		}


		//MPI_Type_free(&dt_recv);
		MPI_Type_free(&dt_send);

		count = 0;

		// problem with the interbuff free !!!
		free(interbuff_offset);
		free(recv_index);
		free(recv_size_by_proc);
		free(send_size_by_proc);

	}
	//free(interbuff_offset);
	free(data_offsets2);
}

void bruck_size(int rank, int num_proc, size_t local_readNum, size_t* number_of_reads_by_procs, int ** data_size, int *new_rank, int *new_size)
{
	MPI_Comm comm = MPI_COMM_WORLD;

	int k, m, j, srank, rrank;
	MPI_Datatype dt_send;
	//MPI_Datatype dt_recv;
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

	int **data_size2 = (int **)malloc(sizeof(int *)*num_proc);
	//we initialize data_offsets
	for(m = 0; m < num_proc; m++){
		data_size[m] = NULL;
		data_size2[m] = calloc( number_of_reads_by_procs[m], sizeof(int));
	}

	size_t *read_by_proc = (size_t *)calloc(num_proc, sizeof(size_t));

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

	//int* interbuff_offset = NULL;
	for(k=1; k<num_proc; k<<=1)
	{
		srank = (rank - k + num_proc) % num_proc;	//Rank to send to
		rrank = (rank + k) % num_proc;	//Rank to recv from

		count = create_send_datatype_for_size(rank, num_proc, number_of_reads_by_procs,
				data_size, k, &dt_send, &recv_index);

		MPI_Type_commit(&dt_send);

		send_size_by_proc = (size_t*)malloc(count*sizeof(size_t));
		recv_size_by_proc = (size_t*)malloc(count*sizeof(size_t));

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

				if (recv_size_by_proc[m] != 0)
					data_size[recv_index[m]] = (int *)malloc(sizeof(int)*(recv_size_by_proc[m]));

				//data_size[recv_index[m]] = realloc(data_size[recv_index[m]], sizeof(int)*(recv_size_by_proc[m]) );
				if (data_size[recv_index[m]])
					data_size[recv_index[m]][0] = 0;

			}
		}

		/*
			MPI_Aint indices_offset[count];
			int blocklens_offset[count];
			MPI_Datatype oldtypes_offset[count];

			for (m = 0; m < count; m++){

				blocklens_offset[m] = (int)(recv_size_by_proc[m]);
				MPI_Get_address(data_offsets[recv_index[m]], &indices_offset[m]);
				oldtypes_offset[m] = MPI_LONG_LONG_INT;
				number_of_reads_by_procs[recv_index[m]] = recv_size_by_proc[m];
			}

			//Create structure of recieve type
			MPI_Type_create_struct(1, blocklens_offset, indices_offset, oldtypes_offset, &dt_recv);
			MPI_Type_commit(&dt_recv);



			int pos=0;
			MPI_Unpack(interbuff_offset, total*8, &pos, MPI_BOTTOM, 1, dt_recv, comm);
		 */

		int *tmp_var = interbuff_offset;

		for (m = 0; m < count; m++){

			memcpy(data_size[recv_index[m]], tmp_var, recv_size_by_proc[m] * sizeof(int));
			tmp_var += recv_size_by_proc[m];
			number_of_reads_by_procs[recv_index[m]] = recv_size_by_proc[m];

		}

		//MPI_Type_free(&dt_recv);
		MPI_Type_free(&dt_send);

		count = 0;

		// problem with the interbuff free !!!
		//if (interbuff_offset)
		free(interbuff_offset);
		free(recv_index);
		free(recv_size_by_proc);
		free(send_size_by_proc);


	}
	free(data_size2);
}

void read_data_for_writing(int rank, int num_proc, size_t local_readNum, char *file_name,
		size_t *number_of_reads_by_procs, size_t *buffs_by_procs, char *** data,
		int *new_rank, int *new_size, size_t *new_offset, MPI_File in, MPI_Info finfo, MPI_Comm COMM_WORLD)
{

	/*
	 * task: IN read_data_for_writing
	 */

	size_t k;
	double t;
	size_t new_data_sz = 0;

	//MPI_File in;
	char **data2 = NULL;

	MPI_Datatype dt_data;
	MPI_Datatype dt_view;

	for (k = 0; k < local_readNum; k++){
		new_data_sz += new_size[k];
	}
	//The data in which what we read will be kept
	data2 = (char**)malloc( num_proc * sizeof(char*));

	// we compute the size of
	// data for each rank and we put it in
	// buffs and buffs_by_proc is
	// the size of the buffer to send
	size_t *buffs_by_procs2 = (size_t*)calloc( num_proc, sizeof(size_t));
	size_t m = 0;

	for(m = 0; m < local_readNum; m++)
	{
		buffs_by_procs2[new_rank[m]] += new_size[m];
		number_of_reads_by_procs[new_rank[m]]++;
	}

	for(m = 0; m < num_proc; m++)
	{
		buffs_by_procs[(rank + m)%num_proc] = buffs_by_procs2[(rank - m + num_proc)%num_proc];
	}

	free(buffs_by_procs2);

	//Allocate data
	for(m = 0; m < num_proc; m++)
	{
		data2[m] = (char*)malloc(buffs_by_procs[m]*sizeof(char) + 1);
		data2[m][buffs_by_procs[m]] = 0;
	}

	//All is here! Creates the datatype for reading
	create_read_dt(rank, num_proc, new_rank, new_size, data2, &dt_data, local_readNum);
	MPI_Type_commit(&dt_data);

	//fprintf(stderr, "%d ::::: [read_data_for_writing] local_readNum = %zu \n", rank, local_readNum);
	//fprintf(stderr, "%d ::::: [read_data_for_writing] data size = %zu \n", rank,  new_data_sz);

	//indexed_data_type_phase2 is for the view and contains the source offset sorted

	MPI_Type_create_hindexed(local_readNum, new_size, (MPI_Aint*)new_offset, MPI_CHAR, &dt_view);
	MPI_Type_commit(&dt_view);

	//TODO: see if initialization is needed
	t = MPI_Wtime();

	MPI_File_set_view(in, 0, MPI_CHAR, dt_view, "native", finfo);
	MPI_File_read(in, MPI_BOTTOM, 1, dt_data, MPI_STATUS_IGNORE);

	//MPI_File_read_at(in, new_offset[0], MPI_BOTTOM, new_data_sz, dt_data, MPI_STATUS_IGNORE);
	//MPI_File_read_all(in, MPI_BOTTOM, 1, dt_data, MPI_STATUS_IGNORE);
	MPI_Barrier(COMM_WORLD);

	//we don't need information of input offset
	MPI_Type_free(&dt_data);
	MPI_Type_free(&dt_view);
	//Close file
	//MPI_File_close(&in);
	/*
	int j =0;
	for (j = 0; j < num_proc; j++){
		fprintf(stderr, "rank %d ::::: [read_data_for_writing] data2 = \n", rank);
		for (m = 0; m < buffs_by_procs[j]; m++){
			fprintf(stderr, "%c", data2[j][m]);
		}
		fprintf(stderr, "\n");
	}
	*/

	*data = data2;
}

void send_size_t_all_to_master(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs, size_t *start_size_per_job, size_t *all_data, size_t *data)
{
	MPI_Status status;
	int j, k;
	if (rank == master){
		//we copy element for rank master_2
		size_t st = start_size_per_job[master];
		for (k = 0; k < size_per_jobs[master]; k++){
			all_data[st] = data[k];
			st++;
		}

		for(j = 0; j < num_proc; j++){

			if (j != master){

				size_t *temp_buf =(size_t *) malloc(size_per_jobs[j]* sizeof(size_t));
				MPI_Recv(temp_buf, size_per_jobs[j], MPI_LONG_LONG_INT, j, 0, MPI_COMM_WORLD, &status);

				size_t st = start_size_per_job[j];
				for (k = 0; k < size_per_jobs[j]; k++){
					all_data[st] = temp_buf[k];
					st++;
				}
				free(temp_buf);
			}
		}
	}
	else{

		MPI_Send(data, size, MPI_LONG_LONG_INT, master,  0, MPI_COMM_WORLD);
	}
}

void send_size_t_all_to_master_bitonic(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs,
		size_t *start_size_per_job, size_t *all_data, size_t *data, size_t start_index)
{
	MPI_Status status;
	int j, k;

	if (rank == master){
		// we copy element for rank master_2
		// we eliminate zeros from the
		// beginning of the vector
		size_t st = start_size_per_job[master];
		for (k = 0; k < size_per_jobs[master]; k++){
			all_data[st] = data[k+start_index];
			st++;
		}

		for(j = 0; j < num_proc; j++){

			if (j != master){

				size_t *temp_buf =(size_t *) malloc(size_per_jobs[j]* sizeof(size_t));
				MPI_Recv(temp_buf, size_per_jobs[j], MPI_LONG_LONG_INT, j, 0, MPI_COMM_WORLD, &status);

				size_t st = start_size_per_job[j];
				for (k = 0; k < size_per_jobs[j]; k++){
					all_data[st] = temp_buf[k];
					st++;
				}
				free(temp_buf);
			}
		}
	}
	else{

		MPI_Send(data, size, MPI_LONG_LONG_INT, master,  0, MPI_COMM_WORLD);
	}
}

void send_size_t_all_to_master_bitonic_V2(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs,
		size_t *start_size_per_job, size_t *all_data, size_t *data, size_t start_index)
{
	MPI_Status status;
	int j, k;

	if (rank == master){
		// we copy element for rank master_2
		// we eliminate zeros from the
		// beginning of the vector

		size_t st = start_size_per_job[master];
		for (k = 0; k < size_per_jobs[master]; k++){
			all_data[st] = data[k+start_index];
			st++;
		}

		for(j = 0; j < num_proc; j++){

			if (j != master){

				size_t *temp_buf =(size_t *) malloc(size_per_jobs[j]* sizeof(size_t));
				MPI_Recv(temp_buf, size_per_jobs[j], MPI_LONG_LONG_INT, j, 0, MPI_COMM_WORLD, &status);
				/*
				fprintf(stderr, "%d ::::: [send_size_t_all_to_master_bitonic_V2] rank %d recieve from %d size = %zu \n",
								rank, rank, j, size_per_jobs[j]);
				*/
				size_t st = start_size_per_job[j];
				/*
				fprintf(stderr, "%d ::::: [send_size_t_all_to_master_bitonic_V2] rank %d copy %zu from %zu in all_data at start %zu \n",
												rank, rank, size_per_jobs[j], j, start_size_per_job[j]);
				*/
				for (k = 0; k < size_per_jobs[j]; k++){
					all_data[st] = temp_buf[k];
					st++;
				}
				free(temp_buf);
			}
		}
	}
	else{
		//fprintf(stderr, "%d ::::: [send_size_t_all_to_master_bitonic_V2] rank %d send to %d at start index= %zu size = %zu \n",rank,
		//		rank, master, start_index, size);

		MPI_Send(data + start_index, size, MPI_LONG_LONG_INT, master,  0, MPI_COMM_WORLD);
	}
}



void send_size_t_master_to_all(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs, size_t *start_size_per_job, size_t *all_data, size_t *data)
{
	MPI_Status status;
	int j, k;
	if (rank != master){
		//fprintf(stderr, "%d ::::: [send_size_t_master_to_all] rank %d recv %zu from %d \n",rank, rank, size, master);

		MPI_Recv(data, size, MPI_LONG_LONG_INT, master, 0, MPI_COMM_WORLD, &status);


	}
	else {

		size_t ind = start_size_per_job[master];
			for (k = 0; k < (size_per_jobs[master]); k++){
			data[k] = all_data[ind];
			ind++;
		}

		for(j = 0; j < num_proc; j++){

			if (j != master){
				//fprintf(stderr, "%d ::::: [send_size_t_master_to_all] rank %d send %zu to %d from %zu\n",
				//	rank, rank, size_per_jobs[j], j, start_size_per_job[j]);
				MPI_Send(&all_data[start_size_per_job[j]],
						size_per_jobs[j], MPI_LONG_LONG_INT, j, 0, MPI_COMM_WORLD);

			}
		}
	}
}

void send_int_all_to_master(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs, size_t *start_size_per_job, int *all_data, int *data)
{
	MPI_Status status;
	int j, k;
	if (rank == master){
		//we copy element for rank master_2
		size_t st = start_size_per_job[master];
		for (k = 0; k < size_per_jobs[master]; k++){
			all_data[st] = data[k];
			st++;
		}

		for(j = 0; j < num_proc; j++){

			if (j != master){

				int *temp_buf =(int *) malloc(size_per_jobs[j]* sizeof(int));
				MPI_Recv(temp_buf, size_per_jobs[j], MPI_INT, j, 0, MPI_COMM_WORLD, &status);

				size_t st = start_size_per_job[j];
				for (k = 0; k < size_per_jobs[j]; k++){
					all_data[st] = temp_buf[k];
					st++;
				}
				free(temp_buf);
			}
		}
	}
	else{

		MPI_Send(data, size, MPI_INT, master,  0, MPI_COMM_WORLD);
	}
}

void send_int_master_to_all(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs, size_t *start_size_per_job, int *all_data, int *data)
{
	MPI_Status status;
	int j, k;

	if (rank != master){
		MPI_Recv(data, size, MPI_INT, master, 0, MPI_COMM_WORLD, &status);
	}
	else {

		size_t ind = start_size_per_job[master];
		for (k = 0; k < size_per_jobs[master]; k++){
			data[k] = all_data[ind];
			ind++;
		}

		for(j = 0; j < num_proc; j++){

			if (j != master){

				MPI_Send(&all_data[start_size_per_job[j]],
						size_per_jobs[j], MPI_INT, j, 0, MPI_COMM_WORLD);
			}
		}
	}
}

void writeSam(int rank, char* output_dir, char* header, size_t local_readNum, char* chrName, Read* chr,
		size_t* offset_dest, int num_proc, MPI_Comm COMM_WORLD, char *file_name, MPI_File in, MPI_Info finfo,
		int compression_level){

	//fprintf(stderr, "rank %d in write sam step 1 local_readNum = %zu \n",rank, local_readNum);
	//task: variables declaration
	size_t j;
	size_t k;
	int ierr;

	MPI_Status status;
	size_t dataSize;

	//size_t *all_offset_dest_index_phase1;
	int *all_rank_to_send=NULL;
	size_t *all_offset_dest_file_to_send=NULL;
	int *all_read_size_to_send=NULL;
	size_t *all_offset_source_file_to_send=NULL;

	/*
	 * phase 1 variables
	 */
	size_t *new_offset_dest_phase1 = (size_t *) malloc(local_readNum* sizeof(size_t));
	int *new_rank_phase1 = (int *) malloc(local_readNum* sizeof(int));
	size_t *new_offset_source_phase1 = (size_t *) malloc(local_readNum* sizeof(size_t));
	//size_t *new_offset_source_index_phase1 = (size_t*)malloc(local_readNum*sizeof(size_t));
	int *new_read_size_phase1 = (int *) malloc(local_readNum* sizeof(int));
	size_t *all_offset_dest_sorted_index_phase1=NULL;
	size_t *pbs_local_dest_offset_index=NULL;
	size_t *pbs_local_dest_offset=NULL;

	/*
	 * phase 2 variables
	 */

	size_t *all_offset_source_file_to_send_phase2=NULL;
	int *all_read_size_phase2_to_send=NULL;
	size_t *all_offset_dest_file_phase2_to_send=NULL;
	int *all_rank_phase2_to_send=NULL;

	//vector allocation for the phase 2
	size_t *new_offset_dest_phase2 = (size_t *) malloc(local_readNum* sizeof(size_t));
	size_t *new_offset_source_sorted_phase2 = (size_t *) malloc(local_readNum* sizeof(size_t));
	int *new_read_size_phase2 = (int *) malloc(local_readNum* sizeof(int));
	int *new_read_size_sorted_phase3 = NULL;// = (int *) malloc(local_readNum* sizeof(int));
	int *new_rank_phase2 = (int *) malloc(local_readNum* sizeof(int));
	size_t *new_offset_dest_index_phase2 = (size_t *) malloc(local_readNum* sizeof(size_t));
	size_t *pbs_local_source_offset=NULL;
	size_t *pbs_local_source_offset_index=NULL;

	size_t *all_offset_source_sorted_index_phase2=NULL;
	//variables for the writing part
	//new_offset_source_index_phase1[0] = 0;

	//the MPI datatype
	MPI_Datatype Datatype_Read_to_write;

	//variables for MPI writes and read
	MPI_File out;
	char* path;
	double time_count;

	size_t *offset_source;
	offset_source = (size_t*)malloc(local_readNum*sizeof(size_t));
	offset_source[0] = 0;

	int *size_source;
	size_source = (int*)malloc(local_readNum*sizeof(int));
	size_source[0] = 0;

	/* TODO Change master_job at each iteration
	 * use master_rank_global_count
	 */
	int master_job_phase_1 = 0;
	int master_job_phase_2 = 0;
	dataSize = 0;

	char **data2;

	//task Init offset and size for source - free chr
	Read *read_tmp = chr;
	dataSize = init_offset_and_size_free_chr(offset_source, size_source, read_tmp, local_readNum);

	/*
	 * before sending all the offsets and size
	 * we sort the offset of destination locally
	 * offset_dest_file_phase1
	 * offset_source_file_phase1
	 * read_size_phase1
	 * task : QSORT PHASE 1
	*/
	size_t *offset_dest_sorted;
	int *size_source_sorted;
	size_t *offset_source_sorted;

	//init indices for qksort
	size_t *offset_dest_index = (size_t*)malloc(local_readNum*sizeof(size_t));

	for(j = 0; j < local_readNum; j++){
		offset_dest_index[j] = j;
	}

	//qksort
	base_arr2 = offset_dest;
	time_count = MPI_Wtime();
	qksort(offset_dest_index, local_readNum, sizeof(size_t), 0, local_readNum - 1, compare_size_t);
	fprintf(stderr, "Rank %d :::::[WRITE] Time to sort local offset source %f seconds\n", rank, MPI_Wtime() - time_count);

	offset_dest_sorted = (size_t*)malloc(local_readNum*sizeof(size_t));
	size_source_sorted = (int*)malloc(local_readNum*sizeof(int));
	offset_source_sorted = (size_t*)malloc(local_readNum*sizeof(size_t));

	//Get result
	for(j = 0; j < local_readNum; j++){
		offset_dest_sorted[j] = offset_dest[offset_dest_index[j]];
		offset_source_sorted[j] = offset_source[offset_dest_index[j]];
		size_source_sorted[j] = size_source[offset_dest_index[j]];

	}

	for ( j = 0; j < local_readNum-1; j++){
				assert( offset_dest_sorted[j] < offset_dest_sorted[j + 1]);
	}

	free(offset_dest);
	free(offset_source);
	free(size_source);
	free(offset_dest_index);

	/* task  Phase 1: Gather all destination offset
	 *
	 *
	 * In this phase we are going to sort the destination
	 * offset in the destination file. In order to
	 * build consecutive blocks to write for each
	 * jobs.
	 *
	 */

	//first we get the vector of all offset in destination file
	size_t *all_offset_dest_file_phase1=NULL;
	size_t *all_offset_source_file_phase1=NULL;
	int *all_read_size_phase1=NULL;

	size_t total_num_read = 0;

	MPI_Reduce(&local_readNum, &total_num_read, 1, MPI_LONG_LONG_INT, MPI_SUM, master_job_phase_1, COMM_WORLD);
	MPI_Reduce(&local_readNum, &total_num_read, 1, MPI_LONG_LONG_INT, MPI_SUM, master_job_phase_2, COMM_WORLD);

	if (rank == master_job_phase_1){
		all_offset_dest_file_phase1 = (size_t *) malloc (total_num_read * sizeof(size_t));
		all_offset_source_file_phase1 = (size_t *) malloc (total_num_read * sizeof(size_t));
		all_read_size_phase1 = (int *) malloc (total_num_read * sizeof(int));
	}

	// vector of number of read per jobs
	size_t *num_reads_per_jobs = (size_t *) malloc(num_proc * sizeof(size_t));
	size_t *start_num_reads_per_jobs_phase1 = (size_t *) malloc((num_proc + 1) * sizeof(size_t));

	// job 0 recieves the number
	// of reads of each rank
	// and put it  in a vector
	MPI_Gather(&local_readNum, 1, MPI_LONG_LONG_INT, &num_reads_per_jobs[rank - master_job_phase_1], 1, MPI_LONG_LONG_INT, master_job_phase_1 , COMM_WORLD);
	MPI_Gather(&local_readNum, 1, MPI_LONG_LONG_INT, &num_reads_per_jobs[rank - master_job_phase_2], 1, MPI_LONG_LONG_INT, master_job_phase_2 , COMM_WORLD);

	if (rank == master_job_phase_1){

		start_num_reads_per_jobs_phase1[0] = 0;

		for (k = 1; k < (num_proc +1); k++){
			start_num_reads_per_jobs_phase1[k] = num_reads_per_jobs[k-1];
		}

		for (k = 1; k < num_proc; k++){
			size_t tmp = start_num_reads_per_jobs_phase1[k - 1];
			size_t tmp2 = start_num_reads_per_jobs_phase1[k];
			start_num_reads_per_jobs_phase1[k] = tmp + tmp2;
		}
	}

	if (rank == master_job_phase_1){

		//we copy the first elements in
		size_t st = start_num_reads_per_jobs_phase1[master_job_phase_1];
		for (k = 0; k < num_reads_per_jobs[master_job_phase_1]; k++){
			all_offset_dest_file_phase1[st] = offset_dest_sorted[k];
			st++;
		}

		for(j = 0; j < num_proc; j++){
			if(j != master_job_phase_1)
			{
				size_t *temp_buf =(size_t *) malloc(num_reads_per_jobs[j]* sizeof(size_t));
				MPI_Recv(temp_buf, num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, COMM_WORLD, &status);
				size_t st = start_num_reads_per_jobs_phase1[j];
				for (k = 0; k < num_reads_per_jobs[j]; k++){
					all_offset_dest_file_phase1[st] = temp_buf[k];
					st++;
				}
				free(temp_buf);
			}
		}
	}
	else{
		MPI_Send(offset_dest_sorted, local_readNum, MPI_LONG_LONG_INT, master_job_phase_1,  0, COMM_WORLD);
	}

	if (rank == master_job_phase_1){

		//we copy the first elements in
		size_t st = start_num_reads_per_jobs_phase1[master_job_phase_1];
		for (k = 0; k < num_reads_per_jobs[master_job_phase_1]; k++){
			all_offset_source_file_phase1[st] = offset_source_sorted[k];
			st++;
		}

		for(j = 0; j < num_proc; j++){
			if(j == master_job_phase_1)
				continue;
			size_t *temp_buf =(size_t *) malloc(num_reads_per_jobs[j]* sizeof(size_t));
			MPI_Recv(temp_buf, num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, COMM_WORLD, &status);

			size_t st = start_num_reads_per_jobs_phase1[j];
			for (k = 0; k < num_reads_per_jobs[j]; k++){
				all_offset_source_file_phase1[st] = temp_buf[k];
				st++;
			}
			free(temp_buf);
		}
	}
	else{
		MPI_Send(offset_source_sorted, local_readNum, MPI_LONG_LONG_INT, master_job_phase_1,  0, COMM_WORLD);
	}


	if (rank == master_job_phase_1){
		//we copy the first elements in
		size_t st = start_num_reads_per_jobs_phase1[master_job_phase_1];
		for (k = 0; k < num_reads_per_jobs[master_job_phase_1]; k++){
			all_read_size_phase1[st] = size_source_sorted[k];
			st++;
		}

		for(j = 0; j < num_proc; j++){
			if(j == master_job_phase_1)
				continue;
			int *temp_buf =(int *) malloc(num_reads_per_jobs[j]* sizeof(int));
			MPI_Recv(temp_buf, num_reads_per_jobs[j], MPI_INT, j, 0, COMM_WORLD, &status);


			size_t st = start_num_reads_per_jobs_phase1[j];
			for (k = 0; k < num_reads_per_jobs[j]; k++){
				all_read_size_phase1[st] = temp_buf[k];
				st++;
			}
			free(temp_buf);
		}
	}
	else{
		MPI_Send(size_source_sorted, local_readNum, MPI_INT, master_job_phase_1,  0, COMM_WORLD);
	}

	/**************************/
	// We free some variable
	/**************************/

	free(size_source_sorted);
	free(offset_dest_sorted);
	free(offset_source_sorted);

	/*******************************************
	 * task: Phase 1: Sorting all offset of destination
	 ******************************************/

	/*
	 *  task Phase 1 : bitonic sort of destination
	 *
	 */

	/*
	 * In this section we implement a parallel Bitonic sort
	 * algorithm.
	 * Input are
	 * all_read_size_phase1,
	 * all_offset_dest_index_phase1
	 * all_offset_dest_file_phase1
	 *
	 */

	int num_processes;
	// Cube Dimension
	MPI_Comm_size(COMM_WORLD, &num_processes);
	// the master rank compute the number of
	// dimension is the number of processors where we
	// perform the bitonic sort
	// int dimensions = (int)(log2(num_processes));
	// find next ( must be greater) power, and go one back

	// we broadcast the total number of reads to each rank
	MPI_Bcast(&total_num_read, 1, MPI_LONG_LONG_INT, master_job_phase_1, COMM_WORLD );

	int dimensions = 1;
	while (dimensions <= num_processes)
		dimensions <<= 1;

	dimensions >>= 1;

	// we compute the length of the vector to recieve
	// the we split all_offset_dest_index_phase1 among the dimension processors
	size_t pbs_num_offsets_to_recieve = total_num_read/dimensions;
	// for the last processor will recieve the rest of the division
	size_t pbs_num_offsets_to_recieve_left  =  total_num_read - pbs_num_offsets_to_recieve*dimensions;


	while((pbs_num_offsets_to_recieve_left * dimensions) > pbs_num_offsets_to_recieve){
		dimensions >>= 1;
		pbs_num_offsets_to_recieve = total_num_read/dimensions;
		pbs_num_offsets_to_recieve_left  =  total_num_read - pbs_num_offsets_to_recieve*dimensions;
	}

	// we compute a vector of size dimensions which contain the number
	// of reads to send
	size_t *pbs_local_num_read_per_job = (size_t *)malloc(dimensions * sizeof(size_t));

	// each job create a vector with with the length
	// of the offset vector for each job in [0, dimension[
	for (j = 0; j < dimensions ; j++){
		// we add pbs_num_offsets_to_recieve_left
		// because we want the vector to have the size size
		pbs_local_num_read_per_job[j] = pbs_num_offsets_to_recieve + (pbs_num_offsets_to_recieve_left);
	}

	// the lastest rank get the reads that left
	//pbs_local_num_read_per_job[dimensions - 1] += pbs_num_offsets_to_recieve_left;
	size_t *pbs_start_num_offset_per_jobs_phase1 = (size_t *) malloc((dimensions + 1) * sizeof(size_t));

	// the master job compute the start index of the element
	// to dispatch
	if (rank == master_job_phase_1){

		pbs_start_num_offset_per_jobs_phase1[0] = 0;

		for (k = 1; k < (dimensions +1); k++){
			pbs_start_num_offset_per_jobs_phase1[k] = pbs_local_num_read_per_job[k-1];
		}
		for (k = 1; k < dimensions; k++){
			size_t tmp = pbs_start_num_offset_per_jobs_phase1[k - 1];
			size_t tmp2 = pbs_start_num_offset_per_jobs_phase1[k];
			// we remove the left over reads
			pbs_start_num_offset_per_jobs_phase1[k] = tmp + tmp2 - pbs_num_offsets_to_recieve_left;
		}
	}

	// the processors master_job_phase_1 send the destination offset
	// to all the rank in [0-dimension]
	if (rank < dimensions){

		// pbs_local_offset_to_sort is a table containing the unsorted
		// destination offset
		pbs_local_dest_offset= (size_t *)malloc(sizeof(size_t) * pbs_local_num_read_per_job[rank]);
		//now the master send

		if ( rank != master_job_phase_1 ){
				MPI_Recv(pbs_local_dest_offset, pbs_local_num_read_per_job[rank], MPI_LONG_LONG_INT,
						master_job_phase_1, 0, MPI_COMM_WORLD, &status);
		}
		else {
			//first we copy the data from the master job
			size_t ind = pbs_start_num_offset_per_jobs_phase1[master_job_phase_1];

			for (k = 0; k < pbs_local_num_read_per_job[master_job_phase_1]; k++){
				pbs_local_dest_offset[k] = all_offset_dest_file_phase1[ind];
				ind++;
			}

			for(j = 0; j < dimensions; j++){
				if (j != master_job_phase_1){
					MPI_Send(&all_offset_dest_file_phase1[pbs_start_num_offset_per_jobs_phase1[j]],
						pbs_local_num_read_per_job[j], MPI_LONG_LONG_INT, j, 0, MPI_COMM_WORLD);
				}
			}
		}
		// we build pbs_local_dest_offset_index
		pbs_local_dest_offset_index = (size_t *)malloc(pbs_local_num_read_per_job[rank]*sizeof(size_t));

		for (j = 0; j < pbs_local_num_read_per_job[rank]; j++){

			if (rank == master_job_phase_1){
				pbs_local_dest_offset_index[j] = j + pbs_local_num_read_per_job[rank]*rank;
			}
			else{
				pbs_local_dest_offset_index[j] = j + pbs_local_num_read_per_job[rank]*rank - (rank*pbs_num_offsets_to_recieve_left);
			}
		}

		// now each rank from [0, dimension[
		// is going to bitonic sort
		// input are:
		// pbs_local_dest_offset
		// pbs_local_dest_offset_index

		time_count = MPI_Wtime();

		ParallelBitonicSort(rank, dimensions, pbs_local_dest_offset, pbs_local_dest_offset_index,
				pbs_local_num_read_per_job[rank], pbs_num_offsets_to_recieve_left);
		fprintf(stderr, "Rank %d :::::[WRITE] Time in the bitonic sort %f seconds\n", rank, MPI_Wtime() - time_count);

		//we compute a new total number of reads
		size_t total_num_read_after_bitonic_sort = 0;
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

		while (pbs_local_dest_offset[start_index] == 0){
			start_index++;
		}

		pbs_local_num_read_per_job[rank] -=  start_index;

		all_offset_dest_file_to_send = (size_t *)malloc(sizeof(size_t) * total_num_read);
		all_offset_dest_sorted_index_phase1 = (size_t *)malloc(sizeof(size_t) * total_num_read);

		if (rank == master_job_phase_1){

			pbs_start_num_offset_per_jobs_phase1[0] = 0;

			for (k = 1; k < (dimensions +1); k++){
				pbs_start_num_offset_per_jobs_phase1[k] = pbs_local_num_read_per_job[k-1];
			}
			for (k = 1; k < dimensions; k++){
				size_t tmp = pbs_start_num_offset_per_jobs_phase1[k - 1];
				size_t tmp2 = pbs_start_num_offset_per_jobs_phase1[k];
				// we remove the left over reads
				pbs_start_num_offset_per_jobs_phase1[k] = tmp + tmp2;
			}
		}

		time_count = MPI_Wtime();
		// we gather the offset destination sorted
		send_size_t_all_to_master_bitonic(rank, dimensions, master_job_phase_1, pbs_local_num_read_per_job[rank],
			pbs_local_num_read_per_job, pbs_start_num_offset_per_jobs_phase1,
				all_offset_dest_file_to_send, pbs_local_dest_offset, start_index);

		// we gather the offset destination sorted index
		send_size_t_all_to_master_bitonic(rank, dimensions, master_job_phase_1, pbs_local_num_read_per_job[rank],
			pbs_local_num_read_per_job, pbs_start_num_offset_per_jobs_phase1,
				all_offset_dest_sorted_index_phase1, pbs_local_dest_offset_index, start_index);


		fprintf(stderr, "Rank %d :::::[WRITE] Time for gathering offset and index %f seconds\n", rank, MPI_Wtime() - time_count);

		if (rank == master_job_phase_1) {



			// we create a new  vector from all_offset_dest_file_to_send
			// and all_offset_dest_sorted_index_phase1 without the zero at the beginning.
			/*
			 * FOR DEBUG
			 *
			for ( j = 0; j < total_num_read; j++){

				assert (all_offset_dest_file_to_send[j] != 0);
			}

			for ( j = 0; j < total_num_read-1; j++){
				assert( all_offset_dest_file_to_send[j] < all_offset_dest_file_to_send[j + 1]);
			}
			*/
			// now we apply the new index to all the
			all_read_size_to_send = (int*)malloc(total_num_read*sizeof(int));
			all_offset_source_file_to_send = (size_t*)malloc(total_num_read*sizeof(size_t));
			all_rank_to_send = (int*)malloc(total_num_read*sizeof(int));

			for(j = 0; j < total_num_read; j++){
				all_offset_source_file_to_send[j] = all_offset_source_file_phase1[all_offset_dest_sorted_index_phase1[j]];
				all_read_size_to_send[j] = all_read_size_phase1[all_offset_dest_sorted_index_phase1[j]];
			}

			/*
			 * FOR DEBUG
			 *
			for ( j = 0; j < total_num_read; j++){
				assert(all_offset_source_file_to_send[j] != 0);

			}
			 */

			// we initialize all_offset_rank_to_send constains
			// the rank of the sorted read
			size_t total = 0;
			for(j = 0; j < num_proc; j++){

				total += num_reads_per_jobs[j];

				for (k = 0; k < num_reads_per_jobs[j]; k++){
					all_rank_to_send[start_num_reads_per_jobs_phase1[j] + k] = j;
				}
			}
			assert(total == total_num_read);
		} // end if (rank == master_job_phase_1)
	} //end if (rank < dimensions)

	// the vector to use in next step are
	// all_rank_to_send
	// all_read_size_to_send

	if ( rank == master_job_phase_1 ){
		free(all_offset_dest_file_phase1);
		free(all_read_size_phase1);
		free(all_offset_source_file_phase1);
	}

	if (rank < dimensions){
		free(all_offset_dest_sorted_index_phase1);
		free(pbs_local_dest_offset);
		free(pbs_local_dest_offset_index);
		free(pbs_local_num_read_per_job);
	}

	/*
	 * task: PROBLEM FREE
	 * SEE WITH VALGRIND
	 */
	//if (pbs_local_num_read_per_job)
	//	free(pbs_local_num_read_per_job);

	if (pbs_start_num_offset_per_jobs_phase1)
		free(pbs_start_num_offset_per_jobs_phase1);

	/******************************/
	// task: Phase 1: dispatch destination offset sorted
	/******************************/

	// from the tables num_reads_per_jobs
	// and start_num_reads_per_jobs
	// for each rank j we send num_reads_per_jobs[j]
	// from start_num_reads_per_jobs[j]

	MPI_Barrier(COMM_WORLD);

	num_reads_per_jobs[rank] = local_readNum;

	assert(local_readNum == num_reads_per_jobs[rank]);

	//Go for the communication
	send_size_t_master_to_all(rank, num_proc, master_job_phase_1, local_readNum, num_reads_per_jobs, start_num_reads_per_jobs_phase1,
			all_offset_dest_file_to_send, new_offset_dest_phase1);

	send_size_t_master_to_all(rank, num_proc, master_job_phase_1, local_readNum, num_reads_per_jobs, start_num_reads_per_jobs_phase1,
			all_offset_source_file_to_send, new_offset_source_phase1);

	send_int_master_to_all(rank, num_proc, master_job_phase_1, local_readNum, num_reads_per_jobs, start_num_reads_per_jobs_phase1,
			all_read_size_to_send, new_read_size_phase1);

	send_int_master_to_all(rank, num_proc, master_job_phase_1, local_readNum, num_reads_per_jobs, start_num_reads_per_jobs_phase1,
			all_rank_to_send, new_rank_phase1);


	if (rank == master_job_phase_1){
		//we free pointers
		free(all_read_size_to_send);
		free(all_offset_source_file_to_send);
		free(all_rank_to_send);
	}

	if (rank < dimensions)
		free(all_offset_dest_file_to_send);

	if (start_num_reads_per_jobs_phase1)
		free(start_num_reads_per_jobs_phase1);


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

	//vector of all offset in destination / source file / size of reads / ranks
	size_t *all_offset_dest_file_phase2=NULL;
	size_t *all_offset_source_file_phase2=NULL;
	int *all_read_size_phase2=NULL;
	int *all_rank_phase2=NULL;

	/*
	 * Phase 2: master_2 gets the total of reads
	 */


	size_t total_num_read_phase2 = 0;
	MPI_Reduce(&local_readNum, &total_num_read_phase2, 1, MPI_LONG_LONG_INT, MPI_SUM, master_job_phase_2, COMM_WORLD);

	if (rank == master_job_phase_2){
		all_offset_dest_file_phase2 = (size_t *) malloc (total_num_read * sizeof(size_t));
		all_offset_source_file_phase2 = (size_t *) malloc (total_num_read * sizeof(size_t));
		all_read_size_phase2 = (int *) malloc (total_num_read * sizeof(int));
		all_rank_phase2 = (int *) malloc (total_num_read * sizeof(int));
	}


	/*
	 * Phase 2: master_2 defines the following two vectors
	 */
	// vector of number of read per jobs
	size_t *num_reads_per_jobs_phase2 = (size_t *) malloc(num_proc* sizeof(size_t));
	// vector of index that contains the cumulative sum of the number of reads
	size_t *start_num_reads_per_jobs_phase2 = (size_t *) malloc((num_proc + 1)*sizeof(size_t));

	/*
	 * Phase 2: master_2 receives all local_readNum and adds it to a local vector
	 */
	MPI_Gather(&local_readNum, 1, MPI_LONG_LONG_INT, &num_reads_per_jobs_phase2[rank - master_job_phase_2],
			1, MPI_LONG_LONG_INT, master_job_phase_2, COMM_WORLD);


	if (rank == master_job_phase_2){

		start_num_reads_per_jobs_phase2[0] = 0;

		for (k = 1; k < (num_proc +1); k++){
			start_num_reads_per_jobs_phase2[k] = num_reads_per_jobs_phase2[k-1];
		}

		for (k = 1; k < num_proc; k++){
			size_t tmp = start_num_reads_per_jobs_phase2[k - 1];
			size_t tmp2 = start_num_reads_per_jobs_phase2[k];
			start_num_reads_per_jobs_phase2[k] = tmp + tmp2;
		}
	}

	/*
	 * task Phase 2: Send from all to master_2
	 */

	//Phase 2: Send local_offset_dest from all to master_2
	send_size_t_all_to_master(rank, num_proc, master_job_phase_2, local_readNum, num_reads_per_jobs_phase2, start_num_reads_per_jobs_phase2,
			all_offset_dest_file_phase2, new_offset_dest_phase1);

	//Phase 2: Send local_offset_source from all to master_2
	send_size_t_all_to_master(rank, num_proc, master_job_phase_2, local_readNum, num_reads_per_jobs_phase2, start_num_reads_per_jobs_phase2,
			all_offset_source_file_phase2, new_offset_source_phase1);

	//Phase 2: Send local_size from all to master_2
	send_int_all_to_master(rank, num_proc, master_job_phase_2, local_readNum, num_reads_per_jobs_phase2, start_num_reads_per_jobs_phase2,
			all_read_size_phase2, new_read_size_phase1);

	//Phase 2: Send local_rank from all to master_2
	send_int_all_to_master(rank, num_proc, master_job_phase_2, local_readNum, num_reads_per_jobs_phase2, start_num_reads_per_jobs_phase2,
			all_rank_phase2, new_rank_phase1);



	/**************************/
	// We free some variable
	/**************************/

	free(new_read_size_phase1);
	free(new_offset_dest_phase1);
	free(new_offset_source_phase1);
	free(new_rank_phase1);

	free(num_reads_per_jobs_phase2);

	/*
	 * task phase 2: bitonic sort of source offset
	 */

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
	// we compute the length of the vector to recieve
	// the we split all_offset_dest_index_phase1 among the dimension processors
	pbs_num_offsets_to_recieve = total_num_read/dimensions;
	// for the last processor will recieve the rest of the division
	pbs_num_offsets_to_recieve_left  =  total_num_read - pbs_num_offsets_to_recieve*dimensions;
	//fprintf(stderr, "pbs_num_offsets_to_recieve_left =  %zu \n", pbs_num_offsets_to_recieve_left);

	// we compute a vector of size dimensions which contain the number
	// of reads to send
	size_t *pbs_local_num_read_per_job_phase2 = (size_t *)malloc(dimensions * sizeof(size_t));

	// each job create a vector with with the length
	// of the offset vector for each job in [0, dimension[
	for (j = 0; j < dimensions ; j++){
		// we add pbs_num_offsets_to_recieve_left
		// because we want the vector to have the same size
		pbs_local_num_read_per_job_phase2[j] = pbs_num_offsets_to_recieve + (pbs_num_offsets_to_recieve_left);
	}

	// the lastest rank get the reads that left
	//pbs_local_num_read_per_job[dimensions - 1] += pbs_num_offsets_to_recieve_left;
	size_t *pbs_start_num_offset_per_jobs_phase2 = (size_t *) malloc((dimensions + 1) * sizeof(size_t));

	while((pbs_num_offsets_to_recieve_left * dimensions) > pbs_num_offsets_to_recieve){
		dimensions >>= 1;
		pbs_num_offsets_to_recieve = total_num_read/dimensions;
		pbs_num_offsets_to_recieve_left  =  total_num_read - pbs_num_offsets_to_recieve*dimensions;
	}

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

	// the processors master_job_phase_2 send the destination offset
	// to all the rank in [0-dimension]
	if (rank < dimensions){

		// pbs_local_offset_to_sort is a table containing the unsorted
		// destination offset
		pbs_local_source_offset = (size_t *)malloc(sizeof(size_t) * pbs_local_num_read_per_job_phase2[rank]);
		//now the master send

		if ( rank != master_job_phase_2 ){
				MPI_Recv(pbs_local_source_offset, pbs_local_num_read_per_job_phase2[rank], MPI_LONG_LONG_INT,
						master_job_phase_2, 0, MPI_COMM_WORLD, &status);

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
					MPI_Send(&all_offset_source_file_phase2[pbs_start_num_offset_per_jobs_phase2[j]],
						pbs_local_num_read_per_job_phase2[j], MPI_LONG_LONG_INT, j, 0, MPI_COMM_WORLD);
				}
			}
		}
		// we build pbs_local_dest_offset_index
		pbs_local_source_offset_index = (size_t *)malloc(pbs_local_num_read_per_job_phase2[rank]*sizeof(size_t));

		//fprintf(stderr, "Rank %d :::::[WRITE] pbs_num_offsets_to_recieve_left = %zu \n", rank, pbs_num_offsets_to_recieve_left);

		for (j = 0; j < pbs_local_num_read_per_job_phase2[rank]; j++){

			if (rank == master_job_phase_2){
				pbs_local_source_offset_index[j] = j + pbs_local_num_read_per_job_phase2[rank]*rank;
			}
			else{
				pbs_local_source_offset_index[j] = j + pbs_local_num_read_per_job_phase2[rank]*rank - (rank*pbs_num_offsets_to_recieve_left);
			}
		}

		// now each rank from [0, dimension[
		// is going to bitonic sort
		// input are:
		// pbs_local_dest_offset
		// pbs_local_dest_offset_index

		// we call the parallel bitonic sort
		time_count = MPI_Wtime();
		ParallelBitonicSort(rank, dimensions, pbs_local_source_offset, pbs_local_source_offset_index,
				pbs_local_num_read_per_job_phase2[rank], pbs_num_offsets_to_recieve_left);

		fprintf(stderr, "Rank %d :::::[WRITE][PHASE 2] Time in the bitonic sort %f seconds\n", rank,
				MPI_Wtime() - time_count);

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

		all_offset_source_file_to_send_phase2 = (size_t *)malloc(sizeof(size_t) * total_num_read_phase2);
		all_offset_source_sorted_index_phase2 = (size_t *)malloc(sizeof(size_t) * total_num_read_phase2);

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
		send_size_t_all_to_master_bitonic_V2(rank, dimensions, master_job_phase_2, pbs_local_num_read_per_job_phase2[rank],
				pbs_local_num_read_per_job_phase2, pbs_start_num_offset_per_jobs_phase2,
					all_offset_source_file_to_send_phase2, pbs_local_source_offset, start_index);

		// we gather the offset source sorted index
		send_size_t_all_to_master_bitonic_V2(rank, dimensions, master_job_phase_2, pbs_local_num_read_per_job_phase2[rank],
				pbs_local_num_read_per_job_phase2, pbs_start_num_offset_per_jobs_phase2,
					all_offset_source_sorted_index_phase2, pbs_local_source_offset_index, start_index);

		fprintf(stderr, "Rank %d :::::[WRITE][PHASE 2] Time for gathering offset and index %f seconds\n", rank,
				MPI_Wtime() - time_count);

		if (rank == master_job_phase_2) {

			// we create a new  vector from all_offset_dest_file_to_send
			// and all_offset_dest_sorted_index_phase1 without the zero at the beginning.
			/*
			 * FOR DEBUG
			 *

			for ( j = 0; j < total_num_read; j++){
				assert ( all_offset_source_file_to_send_phase2[j] != 0 );
			}

			for ( j = 0; j < total_num_read-1; j++){
				assert( all_offset_source_file_to_send_phase2[j] < all_offset_source_file_to_send_phase2[j + 1] );

			}
			*/

			// now we apply the new index to all the
			all_read_size_phase2_to_send = (int*)malloc(total_num_read_phase2*sizeof(int));
			all_offset_dest_file_phase2_to_send = (size_t*)malloc(total_num_read_phase2*sizeof(size_t));
			all_rank_phase2_to_send = (int*)malloc(total_num_read_phase2*sizeof(int));

			for(j = 0; j < total_num_read_phase2; j++){
				all_offset_dest_file_phase2_to_send[j] = all_offset_dest_file_phase2[all_offset_source_sorted_index_phase2[j]];
				all_read_size_phase2_to_send[j] = all_read_size_phase2[all_offset_source_sorted_index_phase2[j]];
				all_rank_phase2_to_send[j] = all_rank_phase2[all_offset_source_sorted_index_phase2[j]];
			}

			/*
			 * FOR DEBUG
			for ( j = 0; j < total_num_read_phase2; j++){
				assert(all_offset_dest_file_phase2_to_send[j] != 0);
			}
			*/

			// we initialize all_offset_rank_to_send constains
			// the rank of the sorted read
			size_t total = 0;
			for(j = 0; j < num_proc; j++){
				total += num_reads_per_jobs[j];
			}
			assert(total == total_num_read_phase2);
		} // end if (rank == master_job_phase_2)

	} //end if (rank < dimensions)

	MPI_Barrier(COMM_WORLD);

	/*
	* FOR DEBUG
	*

	if (rank == master_job_phase_2) {
		for ( j = 0; j < (total_num_read_phase2 - 1); j++){
			assert( all_offset_source_file_to_send_phase2[j] < all_offset_source_file_to_send_phase2[j + 1]);
		}

		for ( j = 0; j < total_num_read_phase2; j++){
			assert(all_offset_dest_file_phase2_to_send[j] != 0);
			assert(all_read_size_phase2_to_send[j] != 0);
			assert(all_rank_phase2_to_send[j] < num_proc);
		}
	}
	*/

	if (rank < dimensions){

		free(pbs_local_source_offset);
		free(pbs_local_source_offset_index);
		free(all_offset_source_sorted_index_phase2);
		free(pbs_local_num_read_per_job_phase2);
		free(pbs_start_num_offset_per_jobs_phase2);
	}

	if (rank == master_job_phase_2){

		free(all_offset_dest_file_phase2);
		free(all_read_size_phase2);
		free(all_rank_phase2);
		free(all_offset_source_file_phase2);
	}

	 // task : PROBLEM FREE
	 //	TO WITH VALGRIND IF THIS FREE ARE CORRECT

	//if (pbs_local_num_read_per_job_phase2)
	//	free(pbs_local_num_read_per_job_phase2);

	//if (pbs_start_num_offset_per_jobs_phase2)
	//	free(pbs_start_num_offset_per_jobs_phase2);

	// the vector to use in next step are
	// all_rank_to_send
	// all_read_size_to_send


	 //  task Phase 2: Dispatch everything

	//Phase 2: Send all_offsets_dest from master_2 to all
	send_size_t_master_to_all(rank, num_proc, master_job_phase_2, local_readNum, num_reads_per_jobs, start_num_reads_per_jobs_phase2,
			all_offset_dest_file_phase2_to_send, new_offset_dest_phase2);

	//Phase 2: Send all_offsets_source from master_2 to all
	send_size_t_master_to_all(rank, num_proc, master_job_phase_2, local_readNum, num_reads_per_jobs, start_num_reads_per_jobs_phase2,
			all_offset_source_file_to_send_phase2, new_offset_source_sorted_phase2);

	//Phase 2: Send all_size from master_2 to all
	send_int_master_to_all(rank, num_proc, master_job_phase_2, local_readNum, num_reads_per_jobs, start_num_reads_per_jobs_phase2,
			all_read_size_phase2_to_send, new_read_size_phase2);

	//Phase 2: Send all_rank from master_2 to all
	send_int_master_to_all(rank, num_proc, master_job_phase_2, local_readNum, num_reads_per_jobs, start_num_reads_per_jobs_phase2,
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

	MPI_Barrier(COMM_WORLD);

	/* *****************************************
	 * task Reading phase
	 * we compute the size of
	 * data for each rank and we put it in
	 * buffs and buffs_by_proc is
	 * the size of the buffer to send
	 *******************************************/
	size_t *number_of_reads_by_procs = (size_t*)calloc( num_proc, sizeof(size_t));
	size_t *buffs_by_procs = (size_t*)calloc( num_proc, sizeof(size_t));

	time_count = MPI_Wtime();
	fprintf(stderr, "Rank %d :::::[WRITE] We call read data for writing \n", rank);

	read_data_for_writing(rank, num_proc, local_readNum,
			file_name, number_of_reads_by_procs, buffs_by_procs,
				&data2, new_rank_phase2, new_read_size_phase2,
					new_offset_source_sorted_phase2, in, finfo, COMM_WORLD);

	free(new_offset_source_sorted_phase2);

	MPI_Barrier(COMM_WORLD);
	fprintf(stderr, "Rank %d :::::[WRITE] Time for reading %f seconds \n", rank, MPI_Wtime() - time_count);


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
	 * 	BEGIN BRUCK PHASE 2    	*
	 *****************************/
	//task Beginning of the Bruck phase
	size_t **data_offsets = (size_t **)malloc(sizeof(size_t *)*num_proc);
	int **data_size = (int **)malloc(sizeof(int *)*num_proc);

	time_count = MPI_Wtime();

	bruckWrite(	rank,
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

	fprintf(stderr, "Rank %d :::::[WRITE] Time for Bruck phases %f seconds\n", rank, MPI_Wtime() - time_count);

	int m;
	local_readNum = 0;

	for (m = 0; m < num_proc; m++){
		local_readNum += number_of_reads_by_procs[m];
	}

	//task Writing phase
	MPI_Aint *reads_address_sorted = (MPI_Aint *)calloc(local_readNum, sizeof(MPI_Aint));

	/*
	 * task: Sort before writing
	 */


	free(new_offset_dest_index_phase2);
	fprintf(stderr, "Rank %d :::::[WRITE] free1\n", rank);

	size_t *new_offset_dest_index_phase3 = (size_t *)malloc(sizeof(size_t) * local_readNum);
	//we sort the offset destination
	new_offset_dest_index_phase3[0] = 0;

	for (k = 1; k < local_readNum; k++){
		new_offset_dest_index_phase3[k] = k;
	}
	//task SORT OFFSET DESTINATION
	//now we sort new_offset_dest_phase2

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
		fprintf(stderr, "Rank %d :::::[WRITE] data_size[%d] = %zu\n", rank, m, data_size[m]);

		if (data_size[m])
			free(data_size[m]);

		j += number_of_reads_by_procs[m];
	}
	fprintf(stderr, "Rank %d :::::[WRITE] free2\n", rank);

	if (data_size)
		free(data_size);
	fprintf(stderr, "Rank %d :::::[WRITE] free3\n", rank);

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
	fprintf(stderr, "Rank %d :::::[WRITE] free4\n", rank);


	if (data_offsets != NULL)
		free(data_offsets);

	fprintf(stderr, "Rank %d :::::[WRITE] free5\n", rank);


	free(number_of_reads_by_procs);
	fprintf(stderr, "Rank %d :::::[WRITE] free6\n", rank);

	base_arr2 = data_offsets_to_sort;
	qksort(new_offset_dest_index_phase3, local_readNum, sizeof(size_t), 0, local_readNum - 1, compare_size_t);

	size_t* offsets_sorted = malloc(sizeof(size_t)*local_readNum);
	for (k = 0; k < local_readNum; k++){
		offsets_sorted[k] = data_offsets_to_sort[new_offset_dest_index_phase3[k]];
	}

	free(data_offsets_to_sort);
	MPI_Datatype *data_table = malloc(local_readNum*sizeof(MPI_Datatype));
	new_read_size_sorted_phase3 = malloc(sizeof(int)*local_readNum);

	/*
	 * data_offsets_to_sort
	 * data_size_to_sort
	 * data_reads_to_sort
	 */

	//for the ring pass phase
	size_t local_size_to_send = 0;

	for(k = 0; k < local_readNum; k++){
		new_read_size_sorted_phase3[k] = data_size_to_sort[new_offset_dest_index_phase3[k]];
		reads_address_sorted[k] = (MPI_Aint)data_reads_to_sort[new_offset_dest_index_phase3[k]];
		local_size_to_send += new_read_size_sorted_phase3[k];
		data_table[k] = MPI_CHAR;
	}


	free(new_offset_dest_index_phase3);
	fprintf(stderr, "Rank %d :::::[WRITE] free7\n", rank);
	MPI_Barrier(COMM_WORLD);

	MPI_Type_create_struct(local_readNum, new_read_size_sorted_phase3, (MPI_Aint*)reads_address_sorted,
			data_table, &Datatype_Read_to_write);

	MPI_Type_commit(&Datatype_Read_to_write);

	/*
	 * In this part and before compress the data
	 * we are going to send data to the next rank
	 * in a ring fashion
	 */

	//Here we are going to send the data to a buffer in the next rank job

	//number of block to send
	int blocksize = 1;

	size_t *y_message_sz = (size_t *)calloc(num_proc,  sizeof(size_t));
	size_t *y_read_num = (size_t *)calloc(num_proc,  sizeof(size_t));

	// copy x into the correct location in y
	y_message_sz[rank * blocksize] = local_size_to_send;
	y_read_num[rank * blocksize] = local_readNum;


	int successor = ( rank + 1 ) % num_proc;
	int predecessor = ( rank - 1 + num_proc ) % num_proc;

	int i=0;
	size_t send_offset;
	size_t recv_offset;

	for (i = 0; i < num_proc - 1 ; i++) {

		send_offset = ( ( rank - i + num_proc ) % num_proc );
		recv_offset = ( ( rank - i - 1 + num_proc ) % num_proc );

		// first we size the size of the buffer to send
		MPI_Send( y_message_sz + send_offset, blocksize , MPI_LONG_LONG_INT, successor, 0, COMM_WORLD);
		MPI_Recv( y_message_sz + recv_offset, blocksize , MPI_LONG_LONG_INT, predecessor, 0, COMM_WORLD, &status);

	}

	//we create a buffer for recieved data
	char *char_buff_uncompressed = malloc(y_message_sz[predecessor] * sizeof(char) + 1);
	char_buff_uncompressed[y_message_sz[predecessor]]='\0';

	//now we send data
	MPI_Sendrecv(MPI_BOTTOM, 1, Datatype_Read_to_write, successor, 0, char_buff_uncompressed, y_message_sz[predecessor],
            MPI_CHAR, predecessor, 0, MPI_COMM_WORLD,  &status);

	fprintf(stderr, "rank %d :::: we recieve from %d a char buffer of size %zu \n", rank, predecessor, strlen(char_buff_uncompressed) );


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

	fprintf(stderr, "Rank %d :::::[WRITE] Time for compressing %f seconds :::: uncompressed size = %d ::: compression size = %d \n",
			rank, MPI_Wtime()-time_count, length, compressed_size);

	free(y_message_sz);

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
		fp_header->cache = kh_init(cache);
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

	MPI_Barrier(COMM_WORLD);

	//we trade the blocks
	//in this phase the predeccor become the succesor
	size_t compressed_sz_to_send = compressed_size;
	size_t compressed_sz_to_recv = 0;

	// copy x into the correct location in y
	int predecessor_back = ( rank + 1 ) % num_proc;
	int successor_back = ( rank - 1 + num_proc ) % num_proc;

	// first we size the size of the buffer to send
	MPI_Send( &compressed_sz_to_send, blocksize , MPI_LONG_LONG_INT, successor_back, 0, MPI_COMM_WORLD);
	MPI_Recv( &compressed_sz_to_recv, blocksize , MPI_LONG_LONG_INT, predecessor_back, 0, MPI_COMM_WORLD, &status);

	//we create a buffer for recieved data
	uint8_t *buff_compressed = malloc(compressed_sz_to_recv * sizeof(uint8_t));
	//char *buff_compressed = malloc(compressed_sz_to_recv * sizeof(char));
	//now we send data
	MPI_Sendrecv(compressed_buff, compressed_sz_to_send, MPI_UNSIGNED_CHAR, successor_back, 0,
			buff_compressed, compressed_sz_to_recv, MPI_UNSIGNED_CHAR, predecessor_back, 0, MPI_COMM_WORLD,  &status);

	free(compressed_buff);
	fprintf(stderr, "rank %d :::: we recieve from %d a compressed buffer of size %zu \n", rank,
			successor, compressed_sz_to_recv );

	MPI_Barrier(COMM_WORLD);
	size_t compSize = compressed_sz_to_recv;


	/*
	 * Now we write results of compression
	 */
	MPI_Barrier(COMM_WORLD);
	size_t write_offset = 0;

	MPI_Offset * y = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
	MPI_Offset * y2 = (MPI_Offset *) calloc(num_proc+1, sizeof(MPI_Offset));

	MPI_Gather(&compSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

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

	//we do a gather in replacement of the the ring pass
	fprintf(stderr, "Proc %d:::: we call MPI_Scatter \n", rank);

	MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	fprintf(stderr, "Proc %d:::: we recieve offset %zu\n", rank, (size_t)write_offset);

	/*
	 * now we compute offset where to write
	 */

	// we create the path where to write for collective write
	path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
	sprintf(path, "%s/%s.bam", output_dir, chrName);

	if(!rank)
		fprintf(stderr, "rank %d :::: Opening the file %s \n", rank, path );
	
		
	//task FINE TUNING FINFO FOR WRITING OPERATIONS
	/*
	MPI_Info_set(finfo,"striping_factor","4");
	MPI_Info_set(finfo,"striping_unit","1610612736"); //1G striping

	MPI_Info_set(finfo,"nb_proc","4");
	MPI_Info_set(finfo,"cb_nodes","4");
	MPI_Info_set(finfo,"cb_block_size","1610612736");
	MPI_Info_set(finfo,"cb_buffer_size","1610612736");
	*/

	ierr = MPI_File_open(COMM_WORLD, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);

	if (ierr) {
		fprintf(stderr, "Rank %d failed to open %s.\nAborting.\n\n", rank, path);
		MPI_Abort(COMM_WORLD, ierr);
		exit(2);
	}
	else{
		if(!rank)fprintf(stderr, "%s.bam successfully opened\n", chrName);
	}

	time_count = MPI_Wtime();

	if (rank == 0 ) {
		fprintf(stderr, "Proc rank %d ::: we write the header \n", rank);
		MPI_File_write(out, compressed_header, compressed_size_header, MPI_BYTE, MPI_STATUS_IGNORE);
	}
	free(compressed_header);

	MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
	MPI_File_write_all(out, buff_compressed, (size_t)compSize, MPI_BYTE, &status);



	//task FINE TUNING FINFO BACK TO READING OPERATIONS
	/*
	MPI_Info_set(finfo,"striping_factor","4");
	MPI_Info_set(finfo,"striping_unit","2684354560"); //1G striping

	MPI_Info_set(finfo,"nb_proc","4");
	MPI_Info_set(finfo,"cb_nodes","4");
	MPI_Info_set(finfo,"cb_block_size","2684354560");
	MPI_Info_set(finfo,"cb_buffer_size","2684354560");
	*/

	free(buff_compressed);

	fprintf(stderr, "Rank %d :::::[WRITE] Time for chromosome %s writing %f seconds\n", rank, chrName, MPI_Wtime()-time_count);

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

	/*
	// we create the path where to write for collective write
	path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
	sprintf(path, "%s/outputSam/%s.sam", output_dir, chrName);

	if(!rank)
		fprintf(stderr, "rank %d :::: Opening the file %s \n", rank, path );

	ierr = MPI_File_open(COMM_WORLD, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, MPI_INFO_NULL, &out);

	if (ierr) {
		fprintf(stderr, "Rank %d failed to open %s.\nAborting.\n\n", rank, path);
		MPI_Abort(COMM_WORLD, ierr);
		exit(2);
	}

	else{
		if(!rank)fprintf(stderr, "%s.sam successfully opened\n", chrName);
	}

	if (rank == 0) {
		MPI_File_write(out, header, strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
	}
	//we create a set view with the offset
		//in destination file
		MPI_Datatype file_view;
		MPI_Type_create_hindexed(local_readNum, &new_read_size_sorted_phase3[0], (MPI_Aint*)offsets_sorted, MPI_CHAR, &file_view);
		MPI_Type_commit(&file_view);

	MPI_File_set_view(out, 0, MPI_CHAR, file_view, "native", MPI_INFO_NULL);
	//task WRITE !!
	time_count = MPI_Wtime();

	MPI_File_write( out, MPI_BOTTOM, 1, Datatype_Read_to_write, &status );

	fprintf(stderr, "Rank %d :::::[WRITE] Time for chromosome %s writing %f seconds\n", rank, chrName, MPI_Wtime()-time_count);

	MPI_Barrier(COMM_WORLD);
	MPI_Type_free(&file_view);
	*/
	//MPI_Type_free(&Datatype_Read_to_write);
	MPI_File_close(&out);
	free(path);
	for(m = 0; m < num_proc; m++)
	{
		free(data2[m]);
	}
	if(y)
		free(y);

	if (y2)
		free(y2);

	free(data2);
	free(data_table);
	free(data_reads_to_sort);
	free(reads_address_sorted);
	free(offsets_sorted);
	free(new_read_size_sorted_phase3);
	free(data_size_to_sort);

}

void writeSam_unmapped(int rank, char* output_dir, char* header, size_t local_readNum, char* chrName, Read* chr,
		int num_proc, MPI_Comm COMM_WORLD, char *file_name, MPI_File in, MPI_Info finfo, int compression_level){

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

	int *read_size = (int *) malloc(local_readNum * sizeof(int));
	int *read_size_sorted = (int *)malloc(local_readNum * sizeof(int));
	offset_source_index[0] = 0;

	//the MPI datatype
	MPI_Datatype indexed_datatype_1;

	//variables for MPI writes and read
	//MPI_File in, out;
	MPI_File out;
	char* path;
	double start, finish, io_time, time_count;

	size_t *offset_source;
	offset_source = (size_t*)malloc(local_readNum*sizeof(size_t));
	offset_source[0] = 0;

	int *size_source;
	size_source = (int*)malloc(local_readNum*sizeof(int));
	size_source[0] = 0;
	dataSize = 0;
	char *data;

	int master_job = 0;
	double start_phase2, finish_phase2;

	//we initialize offset source and size_source
	for(j = 0; j < local_readNum; j++){
		 size_source[j] = 0;
		 offset_source[j] = 0;
	}

	int *all_read_size=NULL;
	size_t *all_offset_source_file=NULL;
	size_t *all_offset_source_index=NULL;
	int *all_read_size_to_send=NULL;

	// vector use in the reading and writing
	// for the phase 1
	size_t *new_offset_dest = (size_t *) malloc(local_readNum* sizeof(size_t));
	int *new_rank = (int *) malloc(local_readNum* sizeof(int));
	size_t *new_offset_source = (size_t *) malloc(local_readNum* sizeof(size_t));
	int *new_read_size = (int *) malloc(local_readNum* sizeof(int));

	//MPI_Info finfo;
	//first we are going to read in the source file
	for(j = 0; j < local_readNum; j++){
			//offset is the read size
			new_offset_source[j] = chr->offset_source_file;
			new_read_size[j] = (int)chr->offset; //read length
			dataSize += chr->offset;
			chr = chr->next;
	}


	/******************************************************************************/
	/* 						Phase one:
	 *
	 * In this phase we are going to sort the source
	 * offset. In order to read consecutive blocks each
	 * jobs.
	 *
	 * Each job has new vector of offset read, of offset write
	 * and of read size :: new_offset_source, new_offset_dest,  new_read_size
	 *
	 * We need a new vector with the rank for sending the reads
	 * after reading.
	 *
	 */
	 /*****************************************************************************/

	start_phase2 = MPI_Wtime();

	//first we get the vector of all offset in destination file
	size_t total_num_read = 0;
	MPI_Reduce(&local_readNum, &total_num_read, 1, MPI_LONG_LONG_INT, MPI_SUM, master_job, COMM_WORLD);
	if (rank == master_job){
		all_offset_source_file = (size_t *) malloc (total_num_read * sizeof(size_t));
		all_read_size = (int *) malloc (total_num_read * sizeof(int));
	}


	// vector of number of read per jobs
	size_t *num_reads_per_jobs = (size_t *) malloc(num_proc* sizeof(size_t));
	// vector of index that contains the cumulative sum of the number of reads
	size_t *start_num_reads_per_jobs = (size_t *) malloc((num_proc + 1)*sizeof(size_t));

	// job 1 recieves the number
	// of reads of each rank
	// and put it  in a vector

	MPI_Gather(&local_readNum, 1, MPI_LONG_LONG_INT, &num_reads_per_jobs[rank - master_job], 1, MPI_LONG_LONG_INT, master_job , COMM_WORLD);

	if (rank == master_job){

		start_num_reads_per_jobs[0] = 0;

		for (k = 1; k < (num_proc +1); k++){
			start_num_reads_per_jobs[k] = num_reads_per_jobs[k-1];
		}
		for (k = 1; k < num_proc; k++){
			size_t tmp = start_num_reads_per_jobs[k - 1];
			size_t tmp2 = start_num_reads_per_jobs[k];
			start_num_reads_per_jobs[k] = tmp + tmp2;
		}
	}

	if (rank == master_job){

		//we copy element for rank 1
		size_t st = start_num_reads_per_jobs[master_job];

		for (k = 0; k < num_reads_per_jobs[master_job]; k++){
			all_offset_source_file[st] = new_offset_source[k];
			st++;
		}

		/*
		 * TODO : See if can replace the loop below with
		 * a simple call like:
		 *
		 * MPI_Recv(&all_offset_dest_file[start_num_reads_per_jobs[j]],
		 *  	num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, MPI_COMM_WORLD, &status);
		 *
		 */

		//fprintf(stderr, "Rank %d ::::: Phase 2 :::: recieve from other job", rank);
		for(j = 0; j < num_proc; j++){

			if ( j != master_job ){

				size_t *temp_buf =(size_t *) malloc(num_reads_per_jobs[j]* sizeof(size_t));
				MPI_Recv(temp_buf, num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, COMM_WORLD, &status);

				size_t st = start_num_reads_per_jobs[j];

				for (k = 0; k < num_reads_per_jobs[j]; k++){
					all_offset_source_file[st] = temp_buf[k];
					st++;
				}
				free(temp_buf);
			}
		}

	}
	else{
		MPI_Send(new_offset_source, local_readNum, MPI_LONG_LONG_INT, master_job,  0, COMM_WORLD);
	}

	//we create vector with all the reads sizes
	//attached to the offset in dest file

	if (rank == master_job){

		//we copy element for rank 1
		size_t st = start_num_reads_per_jobs[master_job];
		for (k = 0; k < num_reads_per_jobs[master_job]; k++){
			all_read_size[st] = new_read_size[k];
			st++;
		}

		/*
		 * TODO : See if can replace the loop below with
		 * a simple call like:
		 *
		 * MPI_Recv(&all_offset_dest_file[start_num_reads_per_jobs[j]],
		 *  	num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, MPI_COMM_WORLD, &status);
		 *
		 */

		for(j = 0; j < num_proc; j++){

			if (j != master_job){

				int *temp_buf =(int *) malloc(num_reads_per_jobs[j]* sizeof(int));
				MPI_Recv(temp_buf, num_reads_per_jobs[j], MPI_INT, j, 0, COMM_WORLD, &status);

				size_t st = start_num_reads_per_jobs[j];
				for (k = 0; k < num_reads_per_jobs[j]; k++){
					all_read_size[st] = temp_buf[k];
					st++;
				}

				free(temp_buf);
			}
		}
	}
	else{
		//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we send new_read_size of size = %zu \n", rank, local_readNum);
		MPI_Send(new_read_size, local_readNum, MPI_INT, master_job,  0, COMM_WORLD);
	}

	/**************************/
	// We free some variable
	/**************************/


	free(new_read_size);
	free(new_offset_dest);
	free(new_offset_source);
	free(new_rank);

	//task: Phase two: Sorting all offset sources

	size_t *all_offset_source_to_send;
	if (rank == master_job) {

		all_offset_source_index = (size_t*)malloc(total_num_read*sizeof(size_t));
		//initialize

		for(j = 0; j < total_num_read; j++){
			all_offset_source_index[j] = j;
			assert(all_offset_source_file[k] != 0);
		}

		//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we start sorting offset in source file\n", rank);
		//previous version
		base_arr2 = all_offset_source_file;
		start = MPI_Wtime();
		qksort(all_offset_source_index, total_num_read, sizeof(size_t), 0, total_num_read - 1, compare_size_t);
		finish = MPI_Wtime();
		io_time = finish - start;
		fprintf(stderr, "Rank %d ::::: Phase 2 :::: for chromosome %s time to sort all the offsets of source file = %f \n", rank, chrName ,io_time);

		all_read_size_to_send = (int*)malloc(total_num_read*sizeof(int));
		all_offset_source_to_send = (size_t*)malloc(total_num_read*sizeof(size_t));

		for(j = 0; j < total_num_read; j++){

			all_offset_source_to_send[j] = all_offset_source_file[all_offset_source_index[j]];
			all_read_size_to_send[j] = all_read_size[all_offset_source_index[j]];
		}

		free(all_offset_source_file);
		free(all_read_size);
		free(all_offset_source_index);
	}

	if (rank == master_job) {
		fprintf(stderr, "Rank %d ::::: After sorting source offset total number of read = %zu \n", rank, total_num_read);
	}

	/************************************/
	// we send pieces of vector sorted
	// according to the sources
	// offset to all jobs
	/************************************/


	// from the tables num_reads_per_jobs
	// and start_num_reads_per_jobs
	// for each rank j we send num_reads_per_jobs[j]
	// from start_num_reads_per_jobs[j]


	//task: Phase two: Dispatch all offset sources
	MPI_Barrier(COMM_WORLD);
	/*
	 * first we send the offsets of the reads in
	 * the destination file
	 */

	/*
	 * second we send the offsets of the reads in
	 * the source file
	 */

	if (rank != master_job){

		//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we recieve offset source \n", rank);
		MPI_Recv(offset_source, local_readNum, MPI_LONG_LONG_INT, master_job, 0, COMM_WORLD, &status);

	}
	else {

		size_t ind = start_num_reads_per_jobs[master_job];
		for (k = 0; k < num_reads_per_jobs[master_job]; k++){
			offset_source[k] = all_offset_source_to_send[ind];
			ind++;
		}

		for(j = 0; j < num_proc; j++){

			if (j != master_job){
				//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we send all_offset_sourceto_send \n", rank);
				MPI_Send(&all_offset_source_to_send[start_num_reads_per_jobs[j]], num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, COMM_WORLD);
			}
		}
	}

	/*
	 * three we send the sizes of the reads in
	 *
	 */

	if (rank != master_job){

		//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we recieve new_read_size \n", rank);
		MPI_Recv(size_source, local_readNum, MPI_INT, master_job, 0, COMM_WORLD, &status);
	}
	else {

		size_t ind = start_num_reads_per_jobs[master_job];
		for (k = 0; k < num_reads_per_jobs[master_job]; k++){
			size_source[k] = all_read_size_to_send[ind];
			ind++;
		}
		for(j = 0; j < num_proc; j++){
			if (j != master_job){
				//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we send all_read_size_to_send \n", rank);
				MPI_Send(&all_read_size_to_send[start_num_reads_per_jobs[j]], num_reads_per_jobs[j], MPI_INT, j, 0, COMM_WORLD);
			}
		}
	}

	fprintf(stderr, "Rank %d ::::: Phase 2 :::: finish dispaching all read size and offset  \n", rank);

	MPI_Barrier(COMM_WORLD);

	if (rank == master_job){
		//we free pointers

		free(all_read_size_to_send);
		free(all_offset_source_to_send);
	}

	free(start_num_reads_per_jobs);
	free(num_reads_per_jobs);
	finish_phase2 = MPI_Wtime();
	io_time = finish_phase2 - start_phase2;


	//fprintf(stderr, "Rank %d ::::: Phase 2 ::::for %s time spend in sending and recieving sorted offset for reading %f seconds\n", rank, chrName, io_time);

	//step 1 :: we sort the vector of input offset
	//we sort the offset_sources
	start = MPI_Wtime();

	for(j = 0; j < local_readNum; j++){
		offset_source_index[j] = j;
	}

	start = MPI_Wtime();
	//previous version of local sort with output of permutation
	base_arr2 = offset_source;
	//new version of the local sort
	qksort(offset_source_index, local_readNum, sizeof(size_t), 0, local_readNum - 1, compare_size_t);
	finish = MPI_Wtime();
	io_time = finish - start;
	fprintf(stderr, "Rank %d ::::: for %s time to local sort the input offset = %f seconds\n", rank, chrName, io_time);

	/*
	 * reorder the offset and the size
	 * of the reads according to sort
	 */

	for(j = 0; j < local_readNum; j++){
		offset_source_sorted[j] = offset_source[offset_source_index[j]];
		read_size_sorted[j] = size_source[offset_source_index[j]];
	}

	size_t new_data_sz = 0;
	for (k = 0; k < local_readNum; k++){
			new_data_sz += read_size_sorted[k];
	}

	MPI_Barrier(COMM_WORLD);

	//assert(new_data_sz == dataSize);
	data = (char *)malloc((new_data_sz + 1));
	data[0]=0;


	/*
	 *
	ierr = MPI_File_open(COMM_WORLD, file_name, MPI_MODE_RDONLY, info1, &in);

	MPI_Barrier(COMM_WORLD);
	if (ierr) {
			fprintf(stderr, "Rank %d failed to open the source file.\nAborting.\n\n", rank);
			MPI_Abort(COMM_WORLD, ierr);
			exit(2);
	}
	 */

	MPI_Type_create_hindexed(local_readNum, &read_size_sorted[0],
			(MPI_Aint*)offset_source_sorted, MPI_CHAR, &indexed_datatype_1);
	MPI_Type_commit(&indexed_datatype_1);

	//we open the file
	//TODO: see if initialization is needed
	//data[new_data_sz] = '\0';
	MPI_File_set_view(in, 0, MPI_CHAR, indexed_datatype_1, "native", finfo);
	start = MPI_Wtime();
	MPI_File_read(in, &data[0], new_data_sz, MPI_CHAR, MPI_STATUS_IGNORE);

	//assert(strlen(data) == dataSize);
	finish = MPI_Wtime();
	io_time = finish - start;
	fprintf(stderr, "Rank %d ::::: read source file = %f seconds\n", rank, io_time);

	MPI_Type_free(&indexed_datatype_1);


	//Here we are going to send the data to a buffer in the next rank job
	//The next job will compress data and send it back to the prvious job

	//number of block to send
	int blocksize = 1;

	size_t *y_message_sz = (size_t *)calloc(num_proc,  sizeof(size_t));
	size_t *y_read_num = (size_t *)calloc(num_proc,  sizeof(size_t));

	// copy x into the correct location in y
	y_message_sz[rank * blocksize] = new_data_sz;
	y_read_num[rank * blocksize] = local_readNum;


	int successor = ( rank + 1 ) % num_proc;
	int predecessor = ( rank - 1 + num_proc ) % num_proc;

	int i=0;
	size_t send_offset;
	size_t recv_offset;

	for (i = 0; i < num_proc - 1 ; i++) {

		send_offset = ( ( rank - i + num_proc ) % num_proc );
		recv_offset = ( ( rank - i - 1 + num_proc ) % num_proc );

		// first we size the size of the buffer to send
		MPI_Send( y_message_sz + send_offset, blocksize , MPI_LONG_LONG_INT, successor, 0, MPI_COMM_WORLD);
		MPI_Recv( y_message_sz + recv_offset, blocksize , MPI_LONG_LONG_INT, predecessor, 0, MPI_COMM_WORLD, &status);

	}

	//we create a buffer for recieved data
	char *char_buff_uncompressed = malloc(y_message_sz[predecessor] * sizeof(char) + 1);

	//now we send data
	MPI_Sendrecv(data, new_data_sz, MPI_CHAR, successor, 0, char_buff_uncompressed, y_message_sz[predecessor],
            MPI_CHAR, predecessor, 0, MPI_COMM_WORLD,  &status);

	fprintf(stderr, "rank %d :::: we recieve from %d a char buffer of size %zu \n", rank, predecessor, strlen(char_buff_uncompressed) );


	BGZF *fp;
	fp = calloc(1, sizeof(BGZF));
	int compress_level = 3;
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
	int compressed_size = 0;

	if (fp->uncompressed_block == NULL)
	   fp->uncompressed_block = malloc(fp->uncompressed_block_size);

	input = (void *)char_buff_uncompressed;
	block_length = fp->uncompressed_block_size;
	bytes_written = 0;
	uint8_t *compressed_buff =  malloc(strlen(char_buff_uncompressed) * sizeof(uint8_t));

	fprintf(stderr, "rank %d :::: start loop compression \n", rank);
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
				memcpy(compressed_buff + compressed_size, fp->compressed_block, block_length);
	        	compressed_size +=block_length;
	        	fp->block_address += block_length;
			}
	    //}
	}

	fprintf(stderr, "Rank %d :::::[WRITE] Time for compressing %f seconds :::: uncompressed size = %d ::: compression size = %d \n",
			rank, MPI_Wtime()-time_count, length, compressed_size);

	free(y_message_sz);

	//we compress the neader


	BGZF *fp_header;
	fp_header = calloc(1, sizeof(BGZF));
	uint8_t *compressed_header = NULL;
	int compressed_size_header = 0;

	if (rank == 0) {
		int compress_level = 3;
		int block_length = MAX_BLOCK_SIZE;
		int bytes_written;
		int length = strlen(header);

		fp_header->open_mode = 'w';
		fp_header->uncompressed_block_size = MAX_BLOCK_SIZE;
		fp_header->uncompressed_block = malloc(MAX_BLOCK_SIZE);
		fp_header->compressed_block_size = MAX_BLOCK_SIZE;
		fp_header->compressed_block = malloc(MAX_BLOCK_SIZE);
		fp_header->cache_size = 0;
		fp_header->cache = kh_init(cache);
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

	MPI_Barrier(COMM_WORLD);

	//we trade the blocks
	//in this phase the predeccor become the succesor
	size_t compressed_sz_to_send = compressed_size;
	size_t compressed_sz_to_recv = 0;

	// copy x into the correct location in y
	int predecessor_back = ( rank + 1 ) % num_proc;
	int successor_back = ( rank - 1 + num_proc ) % num_proc;

	// first we size the size of the buffer to send
	MPI_Send( &compressed_sz_to_send, blocksize , MPI_LONG_LONG_INT, successor_back, 0, MPI_COMM_WORLD);
	MPI_Recv( &compressed_sz_to_recv, blocksize , MPI_LONG_LONG_INT, predecessor_back, 0, MPI_COMM_WORLD, &status);

	//we create a buffer for recieved data
	uint8_t *buff_compressed = malloc(compressed_sz_to_recv * sizeof(uint8_t));
	//char *buff_compressed = malloc(compressed_sz_to_recv * sizeof(char));
	//now we send data
	MPI_Sendrecv(compressed_buff, compressed_sz_to_send, MPI_UNSIGNED_CHAR, successor_back, 0,
			buff_compressed, compressed_sz_to_recv, MPI_UNSIGNED_CHAR, predecessor_back, 0, MPI_COMM_WORLD,  &status);

	//fprintf(stderr, "rank %d :::: we recieve from %d a compressed buffer of size %zu \n", rank,
	//		successor, compressed_sz_to_recv );

	MPI_Barrier(COMM_WORLD);
	size_t compSize = compressed_sz_to_recv;
	/*
	 * Now we write results of compression
	 */

	MPI_Barrier(COMM_WORLD);
	size_t write_offset = 0;

	MPI_Offset * y = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
	MPI_Offset * y2 = (MPI_Offset *) calloc(num_proc+1, sizeof(MPI_Offset));

	MPI_Gather(&compSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

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

	//we do a gather in replacement of the the ring pass
	//fprintf(stderr, "Proc %d:::: we call MPI_Scatter \n", rank);

	MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	//fprintf(stderr, "Proc %d:::: we recieve offset %zu\n", rank, (size_t)write_offset);

	/*
	 * now we compute offset where to write
	 */

	// we create the path where to write for collective write
	path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
	sprintf(path, "%s/%s.bam", output_dir, chrName);

	if(!rank)
		fprintf(stderr, "rank %d :::: Opening the file %s \n", rank, path );

	ierr = MPI_File_open(COMM_WORLD, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);

	if (ierr) {
		fprintf(stderr, "Rank %d failed to open %s.\nAborting.\n\n", rank, path);
		MPI_Abort(COMM_WORLD, ierr);
		exit(2);
	}
	else{
		if(!rank)fprintf(stderr, "%s.bam successfully opened\n", chrName);
	}

	time_count = MPI_Wtime();

	if (rank == 0 ) {
		fprintf(stderr, "Proc rank %d ::: we write the header \n", rank);
		MPI_File_write(out, compressed_header, compressed_size_header, MPI_BYTE, MPI_STATUS_IGNORE);
		//we update write _header
	}
	
	MPI_Barrier(COMM_WORLD);

	//task WRITING OPERATIONS FOR UNMAPPED READS

	MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
	MPI_File_write(out, buff_compressed, (size_t)compSize, MPI_BYTE, &status);

	fprintf(stderr, "Rank %d :::::[WRITE] Time for chromosome %s writing %f seconds\n", rank, chrName, MPI_Wtime()-time_count);

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

	/*
	fprintf(stderr, "Proc %d::::local size to write %zu \n", rank,	dataSize);

	//number of block to send
	MPI_Offset *y = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
	MPI_Offset *y2 = (MPI_Offset *) calloc(num_proc+1, sizeof(MPI_Offset));
	//we wait all processors
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(&dataSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
	//now we make a cumulative sum

	int i1 = 0;

	if (rank ==0){
		for (i1 = 1; i1 < (num_proc + 1); i1++) {
			y2[i1] = y[i1-1];
		}
	}

	i1 = 0;
	if (rank ==0){
		for (i1 = 1; i1 < (num_proc +1); i1++) {
			y2[i1] = y2[i1-1] + y2[i1];
		}
	}

	i1 = 0;
	if (rank ==0){
		for (i1 = 0; i1 < (num_proc +1); i1++) {
			y2[i1] = y2[i1] + strlen(header);
		}
	}

	size_t write_offset = 0;
	MPI_Barrier(MPI_COMM_WORLD);
	//we do a gather in replacement of the the ring pass
	MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	// we create the path where to write for collective write
	path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
	sprintf(path, "%s/outputSam/%s.sam", output_dir, chrName);

	*/

	/*
	 * FOR DEBUG PURPOSE
	 *
	 * MPI_COMM_SELF:: each job write in individual file
	 * ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, info, &out);
	 */
	/*
	ierr = MPI_File_open(MPI_COMM_WORLD, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);

	if (ierr) {
		fprintf(stderr, "Rank %d failed to open %s.\nAborting.\n\n", rank, path);
		MPI_Abort(COMM_WORLD, ierr);
		exit(2);
	}
	else{
		if(!rank)fprintf(stderr, "%s.sam successfully opened\n", chrName);
	}

	start = MPI_Wtime();
	if (rank == 0) {
		MPI_File_write(out, header, strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
	}
	//we create a set view with the offset
	//in destination file
	MPI_File_set_view(out, write_offset, MPI_CHAR, MPI_CHAR, "native", finfo);
	MPI_File_write(out, data, dataSize, MPI_CHAR, MPI_STATUS_IGNORE);
	finish = MPI_Wtime();
	io_time = finish - start;
	fprintf(stderr, "Rank %d ::::: for unmapped reads io_time for writing = %f seconds\n", rank, io_time);
	*/

	MPI_Barrier(COMM_WORLD);
	MPI_File_close(&out);
	//don't close in it will be closed by the discordant part
	//MPI_File_close(&in);
	free(read_size);
	free(read_size_sorted);
	free(offset_source_index);
	free(offset_source_sorted);
	free(offset_source);
	free(size_source);
	free(y);
	free(y2);
	free(path);
	free(data);

}

void writeSam_discordant(int rank, char* output_dir, char* header, size_t local_readNum, char* chrName, Read* chr,
		int num_proc, MPI_Comm COMM_WORLD, char *file_name, MPI_File in, MPI_Info finfo, int compression_level){

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

	int *read_size = (int *) malloc(local_readNum * sizeof(int));
	int *read_size_sorted = (int *)malloc(local_readNum * sizeof(int));
	offset_source_index[0] = 0;

	//the MPI datatype
	MPI_Datatype indexed_datatype_1;

	//variables for MPI writes and read
	//MPI_File in, out;
	MPI_File out;
	char* path;
	double start, finish, io_time, time_count;

	size_t *offset_source;
	offset_source = (size_t*)malloc(local_readNum*sizeof(size_t));
	offset_source[0] = 0;

	int *size_source;
	size_source = (int*)malloc(local_readNum*sizeof(int));
	size_source[0] = 0;
	dataSize = 0;
	char *data;

	int master_job = 0;
	double start_phase2, finish_phase2;

	//we initialize offset source and size_source
	for(j = 0; j < local_readNum; j++){
		 size_source[j] = 0;
		 offset_source[j] = 0;
	}

	int *all_read_size=NULL;
	size_t *all_offset_source_file=NULL;
	size_t *all_offset_source_index=NULL;
	int *all_read_size_to_send=NULL;

	// vector use in the reading and writing
	// for the phase 1
	size_t *new_offset_dest = (size_t *) malloc(local_readNum* sizeof(size_t));
	int *new_rank = (int *) malloc(local_readNum* sizeof(int));
	size_t *new_offset_source = (size_t *) malloc(local_readNum* sizeof(size_t));
	int *new_read_size = (int *) malloc(local_readNum* sizeof(int));

	//MPI_Info finfo;
	//first we are going to read in the source file
	for(j = 0; j < local_readNum; j++){
			//offset is the read size
			new_offset_source[j] = chr->offset_source_file;
			new_read_size[j] = (int)chr->offset; //read length
			dataSize += chr->offset;
			chr = chr->next;
	}


	/******************************************************************************/
	/* 						Phase one:
	 *
	 * In this phase we are going to sort the source
	 * offset. In order to read consecutive blocks each
	 * jobs.
	 *
	 * Each job has new vector of offset read, of offset write
	 * and of read size :: new_offset_source, new_offset_dest,  new_read_size
	 *
	 * We need a new vector with the rank for sending the reads
	 * after reading.
	 *
	 */
	 /*****************************************************************************/

	start_phase2 = MPI_Wtime();

	//first we get the vector of all offset in destination file
	size_t total_num_read = 0;
	MPI_Reduce(&local_readNum, &total_num_read, 1, MPI_LONG_LONG_INT, MPI_SUM, master_job, COMM_WORLD);
	if (rank == master_job){
		all_offset_source_file = (size_t *) malloc (total_num_read * sizeof(size_t));
		all_read_size = (int *) malloc (total_num_read * sizeof(int));
	}


	// vector of number of read per jobs
	size_t *num_reads_per_jobs = (size_t *) malloc(num_proc* sizeof(size_t));
	// vector of index that contains the cumulative sum of the number of reads
	size_t *start_num_reads_per_jobs = (size_t *) malloc((num_proc + 1)*sizeof(size_t));

	// job 1 recieves the number
	// of reads of each rank
	// and put it  in a vector

	MPI_Gather(&local_readNum, 1, MPI_LONG_LONG_INT, &num_reads_per_jobs[rank - master_job], 1, MPI_LONG_LONG_INT, master_job , COMM_WORLD);

	if (rank == master_job){

		start_num_reads_per_jobs[0] = 0;

		for (k = 1; k < (num_proc +1); k++){
			start_num_reads_per_jobs[k] = num_reads_per_jobs[k-1];
		}
		for (k = 1; k < num_proc; k++){
			size_t tmp = start_num_reads_per_jobs[k - 1];
			size_t tmp2 = start_num_reads_per_jobs[k];
			start_num_reads_per_jobs[k] = tmp + tmp2;
		}
	}

	if (rank == master_job){

		//we copy element for rank 1
		size_t st = start_num_reads_per_jobs[master_job];

		for (k = 0; k < num_reads_per_jobs[master_job]; k++){
			all_offset_source_file[st] = new_offset_source[k];
			st++;
		}

		/*
		 * TODO : See if can replace the loop below with
		 * a simple call like:
		 *
		 * MPI_Recv(&all_offset_dest_file[start_num_reads_per_jobs[j]],
		 *  	num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, MPI_COMM_WORLD, &status);
		 *
		 */

		//fprintf(stderr, "Rank %d ::::: Phase 2 :::: recieve from other job", rank);
		for(j = 0; j < num_proc; j++){

			if ( j != master_job ){

				size_t *temp_buf =(size_t *) malloc(num_reads_per_jobs[j]* sizeof(size_t));
				MPI_Recv(temp_buf, num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, COMM_WORLD, &status);

				size_t st = start_num_reads_per_jobs[j];

				for (k = 0; k < num_reads_per_jobs[j]; k++){
					all_offset_source_file[st] = temp_buf[k];
					st++;
				}
				free(temp_buf);
			}
		}

	}
	else{
		MPI_Send(new_offset_source, local_readNum, MPI_LONG_LONG_INT, master_job,  0, COMM_WORLD);
	}

	//we create vector with all the reads sizes
	//attached to the offset in dest file

	if (rank == master_job){

		//we copy element for rank 1
		size_t st = start_num_reads_per_jobs[master_job];
		for (k = 0; k < num_reads_per_jobs[master_job]; k++){
			all_read_size[st] = new_read_size[k];
			st++;
		}

		/*
		 * TODO : See if can replace the loop below with
		 * a simple call like:
		 *
		 * MPI_Recv(&all_offset_dest_file[start_num_reads_per_jobs[j]],
		 *  	num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, MPI_COMM_WORLD, &status);
		 *
		 */

		for(j = 0; j < num_proc; j++){

			if (j != master_job){

				int *temp_buf =(int *) malloc(num_reads_per_jobs[j]* sizeof(int));
				MPI_Recv(temp_buf, num_reads_per_jobs[j], MPI_INT, j, 0, COMM_WORLD, &status);

				size_t st = start_num_reads_per_jobs[j];
				for (k = 0; k < num_reads_per_jobs[j]; k++){
					all_read_size[st] = temp_buf[k];
					st++;
				}

				free(temp_buf);
			}
		}
	}
	else{
		//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we send new_read_size of size = %zu \n", rank, local_readNum);
		MPI_Send(new_read_size, local_readNum, MPI_INT, master_job,  0, COMM_WORLD);
	}

	/**************************/
	// We free some variable
	/**************************/


	free(new_read_size);
	free(new_offset_dest);
	free(new_offset_source);
	free(new_rank);

	//task: Phase two: Sorting all offset sources

	size_t *all_offset_source_to_send;
	if (rank == master_job) {

		all_offset_source_index = (size_t*)malloc(total_num_read*sizeof(size_t));
		//initialize

		for(j = 0; j < total_num_read; j++){
			all_offset_source_index[j] = j;
			assert(all_offset_source_file[k] != 0);
		}

		//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we start sorting offset in source file\n", rank);
		//previous version
		base_arr2 = all_offset_source_file;
		start = MPI_Wtime();
		qksort(all_offset_source_index, total_num_read, sizeof(size_t), 0, total_num_read - 1, compare_size_t);
		finish = MPI_Wtime();
		io_time = finish - start;
		fprintf(stderr, "Rank %d ::::: Phase 2 :::: for chromosome %s time to sort all the offsets of source file = %f \n", rank, chrName ,io_time);

		all_read_size_to_send = (int*)malloc(total_num_read*sizeof(int));
		all_offset_source_to_send = (size_t*)malloc(total_num_read*sizeof(size_t));

		for(j = 0; j < total_num_read; j++){

			all_offset_source_to_send[j] = all_offset_source_file[all_offset_source_index[j]];
			all_read_size_to_send[j] = all_read_size[all_offset_source_index[j]];
		}

		free(all_offset_source_file);
		free(all_read_size);
		free(all_offset_source_index);
	}

	if (rank == master_job) {
		fprintf(stderr, "Rank %d ::::: After sorting source offset total number of read = %zu \n", rank, total_num_read);
	}

	/************************************/
	// we send pieces of vector sorted
	// according to the sources
	// offset to all jobs
	/************************************/


	// from the tables num_reads_per_jobs
	// and start_num_reads_per_jobs
	// for each rank j we send num_reads_per_jobs[j]
	// from start_num_reads_per_jobs[j]


	//task: Phase two: Dispatch all offset sources
	MPI_Barrier(COMM_WORLD);
	/*
	 * first we send the offsets of the reads in
	 * the destination file
	 */

	/*
	 * second we send the offsets of the reads in
	 * the source file
	 */

	if (rank != master_job){

		//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we recieve offset source \n", rank);
		MPI_Recv(offset_source, local_readNum, MPI_LONG_LONG_INT, master_job, 0, COMM_WORLD, &status);

	}
	else {

		size_t ind = start_num_reads_per_jobs[master_job];
		for (k = 0; k < num_reads_per_jobs[master_job]; k++){
			offset_source[k] = all_offset_source_to_send[ind];
			ind++;
		}

		for(j = 0; j < num_proc; j++){

			if (j != master_job){
				//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we send all_offset_sourceto_send \n", rank);
				MPI_Send(&all_offset_source_to_send[start_num_reads_per_jobs[j]], num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, COMM_WORLD);
			}
		}
	}

	/*
	 * three we send the sizes of the reads in
	 *
	 */

	if (rank != master_job){

		//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we recieve new_read_size \n", rank);
		MPI_Recv(size_source, local_readNum, MPI_INT, master_job, 0, COMM_WORLD, &status);
	}
	else {

		size_t ind = start_num_reads_per_jobs[master_job];
		for (k = 0; k < num_reads_per_jobs[master_job]; k++){
			size_source[k] = all_read_size_to_send[ind];
			ind++;
		}
		for(j = 0; j < num_proc; j++){
			if (j != master_job){
				//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we send all_read_size_to_send \n", rank);
				MPI_Send(&all_read_size_to_send[start_num_reads_per_jobs[j]], num_reads_per_jobs[j], MPI_INT, j, 0, COMM_WORLD);
			}
		}
	}

	fprintf(stderr, "Rank %d ::::: Phase 2 :::: finish dispaching all read size and offset  \n", rank);

	MPI_Barrier(COMM_WORLD);

	if (rank == master_job){
		//we free pointers

		free(all_read_size_to_send);
		free(all_offset_source_to_send);
	}

	free(start_num_reads_per_jobs);
	free(num_reads_per_jobs);
	finish_phase2 = MPI_Wtime();
	io_time = finish_phase2 - start_phase2;


	//fprintf(stderr, "Rank %d ::::: Phase 2 ::::for %s time spend in sending and recieving sorted offset for reading %f seconds\n", rank, chrName, io_time);

	//step 1 :: we sort the vector of input offset
	//we sort the offset_sources
	start = MPI_Wtime();

	for(j = 0; j < local_readNum; j++){
		offset_source_index[j] = j;
	}

	start = MPI_Wtime();
	//previous version of local sort with output of permutation
	base_arr2 = offset_source;
	//new version of the local sort
	qksort(offset_source_index, local_readNum, sizeof(size_t), 0, local_readNum - 1, compare_size_t);
	finish = MPI_Wtime();
	io_time = finish - start;
	fprintf(stderr, "Rank %d ::::: for %s time to local sort the input offset = %f seconds\n", rank, chrName, io_time);

	/*
	 * reorder the offset and the size
	 * of the reads according to sort
	 */

	for(j = 0; j < local_readNum; j++){
		offset_source_sorted[j] = offset_source[offset_source_index[j]];
		read_size_sorted[j] = size_source[offset_source_index[j]];
	}

	size_t new_data_sz = 0;
	for (k = 0; k < local_readNum; k++){
			new_data_sz += read_size_sorted[k];
	}

	MPI_Barrier(COMM_WORLD);

	//assert(new_data_sz == dataSize);
	data = (char *)malloc((new_data_sz + 1));
	data[0]=0;


	/*
	 *
	ierr = MPI_File_open(COMM_WORLD, file_name, MPI_MODE_RDONLY, info1, &in);

	MPI_Barrier(COMM_WORLD);
	if (ierr) {
			fprintf(stderr, "Rank %d failed to open the source file.\nAborting.\n\n", rank);
			MPI_Abort(COMM_WORLD, ierr);
			exit(2);
	}
	 */

	MPI_Type_create_hindexed(local_readNum, &read_size_sorted[0],
			(MPI_Aint*)offset_source_sorted, MPI_CHAR, &indexed_datatype_1);
	MPI_Type_commit(&indexed_datatype_1);

	//we open the file
	//TODO: see if initialization is needed
	//data[new_data_sz] = '\0';
	MPI_File_set_view(in, 0, MPI_CHAR, indexed_datatype_1, "native", finfo);
	start = MPI_Wtime();
	MPI_File_read(in, &data[0], new_data_sz, MPI_CHAR, MPI_STATUS_IGNORE);

	//assert(strlen(data) == dataSize);
	finish = MPI_Wtime();
	io_time = finish - start;
	fprintf(stderr, "Rank %d ::::: read source file = %f seconds\n", rank, io_time);

	MPI_Type_free(&indexed_datatype_1);


	//Here we are going to send the data to a buffer in the next rank job
	//The next job will compress data and send it back to the prvious job

	//number of block to send
	int blocksize = 1;

	size_t *y_message_sz = (size_t *)calloc(num_proc,  sizeof(size_t));
	size_t *y_read_num = (size_t *)calloc(num_proc,  sizeof(size_t));

	// copy x into the correct location in y
	y_message_sz[rank * blocksize] = new_data_sz;
	y_read_num[rank * blocksize] = local_readNum;


	int successor = ( rank + 1 ) % num_proc;
	int predecessor = ( rank - 1 + num_proc ) % num_proc;

	int i=0;
	size_t send_offset;
	size_t recv_offset;

	for (i = 0; i < num_proc - 1 ; i++) {

		send_offset = ( ( rank - i + num_proc ) % num_proc );
		recv_offset = ( ( rank - i - 1 + num_proc ) % num_proc );

		// first we size the size of the buffer to send
		MPI_Send( y_message_sz + send_offset, blocksize , MPI_LONG_LONG_INT, successor, 0, MPI_COMM_WORLD);
		MPI_Recv( y_message_sz + recv_offset, blocksize , MPI_LONG_LONG_INT, predecessor, 0, MPI_COMM_WORLD, &status);

	}

	//we create a buffer for recieved data
	char *char_buff_uncompressed = malloc(y_message_sz[predecessor] * sizeof(char) + 1);

	//now we send data
	MPI_Sendrecv(data, new_data_sz, MPI_CHAR, successor, 0, char_buff_uncompressed, y_message_sz[predecessor],
            MPI_CHAR, predecessor, 0, MPI_COMM_WORLD,  &status);

	fprintf(stderr, "rank %d :::: we recieve from %d a char buffer of size %zu \n", rank, predecessor, strlen(char_buff_uncompressed) );


	BGZF *fp;
	fp = calloc(1, sizeof(BGZF));
	int compress_level = 3;
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
	int compressed_size = 0;

	if (fp->uncompressed_block == NULL)
	   fp->uncompressed_block = malloc(fp->uncompressed_block_size);

	input = (void *)char_buff_uncompressed;
	block_length = fp->uncompressed_block_size;
	bytes_written = 0;
	uint8_t *compressed_buff =  malloc(strlen(char_buff_uncompressed) * sizeof(uint8_t));

	fprintf(stderr, "rank %d :::: start loop compression \n", rank);
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
				memcpy(compressed_buff + compressed_size, fp->compressed_block, block_length);
	        	compressed_size +=block_length;
	        	fp->block_address += block_length;
			}
	    //}
	}

	fprintf(stderr, "Rank %d :::::[WRITE] Time for compressing %f seconds :::: uncompressed size = %d ::: compression size = %d \n",
			rank, MPI_Wtime()-time_count, length, compressed_size);

	free(y_message_sz);

	//we compress the neader


	BGZF *fp_header;
	fp_header = calloc(1, sizeof(BGZF));
	uint8_t *compressed_header = NULL;
	int compressed_size_header = 0;

	if (rank == 0) {
		int compress_level = 3;
		int block_length = MAX_BLOCK_SIZE;
		int bytes_written;
		int length = strlen(header);

		fp_header->open_mode = 'w';
		fp_header->uncompressed_block_size = MAX_BLOCK_SIZE;
		fp_header->uncompressed_block = malloc(MAX_BLOCK_SIZE);
		fp_header->compressed_block_size = MAX_BLOCK_SIZE;
		fp_header->compressed_block = malloc(MAX_BLOCK_SIZE);
		fp_header->cache_size = 0;
		fp_header->cache = kh_init(cache);
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

	MPI_Barrier(COMM_WORLD);

	//we trade the blocks
	//in this phase the predeccor become the succesor
	size_t compressed_sz_to_send = compressed_size;
	size_t compressed_sz_to_recv = 0;

	// copy x into the correct location in y
	int predecessor_back = ( rank + 1 ) % num_proc;
	int successor_back = ( rank - 1 + num_proc ) % num_proc;

	// first we size the size of the buffer to send
	MPI_Send( &compressed_sz_to_send, blocksize , MPI_LONG_LONG_INT, successor_back, 0, MPI_COMM_WORLD);
	MPI_Recv( &compressed_sz_to_recv, blocksize , MPI_LONG_LONG_INT, predecessor_back, 0, MPI_COMM_WORLD, &status);

	//we create a buffer for recieved data
	uint8_t *buff_compressed = malloc(compressed_sz_to_recv * sizeof(uint8_t));
	//char *buff_compressed = malloc(compressed_sz_to_recv * sizeof(char));
	//now we send data
	MPI_Sendrecv(compressed_buff, compressed_sz_to_send, MPI_UNSIGNED_CHAR, successor_back, 0,
			buff_compressed, compressed_sz_to_recv, MPI_UNSIGNED_CHAR, predecessor_back, 0, MPI_COMM_WORLD,  &status);

	//fprintf(stderr, "rank %d :::: we recieve from %d a compressed buffer of size %zu \n", rank,
	//		successor, compressed_sz_to_recv );

	MPI_Barrier(COMM_WORLD);
	size_t compSize = compressed_sz_to_recv;
	/*
	 * Now we write results of compression
	 */

	MPI_Barrier(COMM_WORLD);
	size_t write_offset = 0;

	MPI_Offset * y = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
	MPI_Offset * y2 = (MPI_Offset *) calloc(num_proc+1, sizeof(MPI_Offset));

	MPI_Gather(&compSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

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

	//we do a gather in replacement of the the ring pass
	//fprintf(stderr, "Proc %d:::: we call MPI_Scatter \n", rank);

	MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	//fprintf(stderr, "Proc %d:::: we recieve offset %zu\n", rank, (size_t)write_offset);

	/*
	 * now we compute offset where to write
	 */

	// we create the path where to write for collective write
	path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
	sprintf(path, "%s/%s.bam", output_dir, chrName);

	if(!rank)
		fprintf(stderr, "rank %d :::: Opening the file %s \n", rank, path );

	ierr = MPI_File_open(COMM_WORLD, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);

	if (ierr) {
		fprintf(stderr, "Rank %d failed to open %s.\nAborting.\n\n", rank, path);
		MPI_Abort(COMM_WORLD, ierr);
		exit(2);
	}
	else{
		if(!rank)fprintf(stderr, "%s.bam successfully opened\n", chrName);
	}

	time_count = MPI_Wtime();

	if (rank == 0 ) {
		fprintf(stderr, "Proc rank %d ::: we write the header \n", rank);
		MPI_File_write(out, compressed_header, compressed_size_header, MPI_BYTE, MPI_STATUS_IGNORE);
		//we update write _header
	}

	MPI_Barrier(COMM_WORLD);

	//task WRITING OPERATIONS FOR UNMAPPED READS

	MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
	MPI_File_write(out, buff_compressed, (size_t)compSize, MPI_BYTE, &status);

	fprintf(stderr, "Rank %d :::::[WRITE] Time for chromosome %s writing %f seconds\n", rank, chrName, MPI_Wtime()-time_count);

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

	/*
	fprintf(stderr, "Proc %d::::local size to write %zu \n", rank,	dataSize);

	//number of block to send
	MPI_Offset *y = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
	MPI_Offset *y2 = (MPI_Offset *) calloc(num_proc+1, sizeof(MPI_Offset));
	//we wait all processors
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(&dataSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
	//now we make a cumulative sum

	int i1 = 0;

	if (rank ==0){
		for (i1 = 1; i1 < (num_proc + 1); i1++) {
			y2[i1] = y[i1-1];
		}
	}

	i1 = 0;
	if (rank ==0){
		for (i1 = 1; i1 < (num_proc +1); i1++) {
			y2[i1] = y2[i1-1] + y2[i1];
		}
	}

	i1 = 0;
	if (rank ==0){
		for (i1 = 0; i1 < (num_proc +1); i1++) {
			y2[i1] = y2[i1] + strlen(header);
		}
	}

	size_t write_offset = 0;
	MPI_Barrier(MPI_COMM_WORLD);
	//we do a gather in replacement of the the ring pass
	MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	// we create the path where to write for collective write
	path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
	sprintf(path, "%s/outputSam/%s.sam", output_dir, chrName);

	*/

	/*
	 * FOR DEBUG PURPOSE
	 *
	 * MPI_COMM_SELF:: each job write in individual file
	 * ierr = MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, info, &out);
	 */
	/*
	ierr = MPI_File_open(MPI_COMM_WORLD, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);

	if (ierr) {
		fprintf(stderr, "Rank %d failed to open %s.\nAborting.\n\n", rank, path);
		MPI_Abort(COMM_WORLD, ierr);
		exit(2);
	}
	else{
		if(!rank)fprintf(stderr, "%s.sam successfully opened\n", chrName);
	}

	start = MPI_Wtime();
	if (rank == 0) {
		MPI_File_write(out, header, strlen(header), MPI_CHAR, MPI_STATUS_IGNORE);
	}
	//we create a set view with the offset
	//in destination file
	MPI_File_set_view(out, write_offset, MPI_CHAR, MPI_CHAR, "native", finfo);
	MPI_File_write(out, data, dataSize, MPI_CHAR, MPI_STATUS_IGNORE);
	finish = MPI_Wtime();
	io_time = finish - start;
	fprintf(stderr, "Rank %d ::::: for unmapped reads io_time for writing = %f seconds\n", rank, io_time);
	*/

	MPI_Barrier(COMM_WORLD);
	MPI_File_close(&out);
	MPI_File_close(&in);
	free(read_size);
	free(read_size_sorted);
	free(offset_source_index);
	free(offset_source_sorted);
	free(offset_source);
	free(size_source);
	free(y);
	free(y2);
	free(path);
	free(data);

}


static int partition(void *data, size_t esize, size_t i, size_t k, int (*compare)(const void *key1, const void *key2)){

	size_t *a = data;
	size_t *pval, *temp;
	size_t r[3];
	/*
	 * allocate value for the partition value and swapping
	 */
	if ((pval = (size_t *)malloc(sizeof(size_t))) == NULL)
		return -1;
	if ((temp = (size_t *)malloc(sizeof(size_t))) == NULL){
		free(pval);
		return -1;
	}

	/*
	 * Use the median-of-three method to find partition value
	 */

	int part = (k - i + 1);
	if(part <= 0)
		part = 1;

	r[0] = (rand() % part +i);
	r[1] = (rand() % part +i);
	r[2] = (rand() % part +i);


	/*
	 * TODO: replace the qsort with issort
	 *
	 * issort(r, 3, sizeof(size_t), compare_size_t_V2);
	 */
	qsort(r, 3, sizeof(size_t),compare_size_t_V2);
	memcpy(pval, &a[r[1]], sizeof(size_t));

	/*
	 * Create 2 partitions around the partition value
	 */
	i--;
	k++;
	while(1) {

		/*
		 * move left until an element is found in the wrong partition
		 */

		do {
			k--;
		} while (compare(&a[k], pval) > 0);

		/*
		 * move right until an element is found in the wrong partition
		 */

		do {
			i++;
		} while (compare(&a[i], pval) < 0);

		if (i >= k){
			/*
			 * break when left and right counter cross
			 */
			break;
		}

		else{
			// swap element under the left and right counters
			memcpy(temp, &a[i], sizeof(size_t));
			memcpy(&a[i], &a[k], sizeof(size_t));
			memcpy(&a[k], temp, sizeof(size_t));
		}
	}

	/*
	 * free the storage allocated for partitioning
	 */
	free(pval);
	free(temp);

	/*
	 * return position dividing the two partition
	 */
	return k;

}

/*
 * -------------------------            qksort          ------------------------------
 */

int qksort(void *data, size_t size, size_t esize, size_t i, size_t k, int (*compare)(const void *key1, const void *key2)){

	size_t j;

	/*
	 * stop recursion when it is not possible to partition further
	 * when calling qksort:
	 * i = 0
	 * k = size-1
	 */

	while (i < k){

		/*
		 * find where to partition the elements
		 */


		if ((j = partition(data, esize, i, k, compare)) < 0){
			return -1;
		}
		/*
		 * recursively sort the left partition
		 */
		if (qksort(data, size, esize, i, j, compare) < 0)
			return -1;
		/*
		 * iterate and sort the right partition
		 */
		i = j + 1;
	}
	return 0;
}


int issort(void *data, size_t size, size_t esize, int (*compare)(const void *key, const void *key2)){

	size_t *a = data;
	size_t *key;
	size_t i,j;

	if ((key = (size_t *)malloc(sizeof(size_t))) == NULL)
		return -1;

	for ( j = 1; j < size; j++){
		memcpy(key, &a[j], sizeof(size_t));
		i = j - 1;
		while (i >= 0 && compare(&a[i], key) > 0){
			memcpy(&a[(i + 1)], &a[i], sizeof(size_t));
			i--;
		}
		memcpy(&a[(i + 1)], key, sizeof(size_t));
	}

	free(key);
	return 0;

}

void create_DT(int rank, int numprocs, int k, MPI_Datatype oldtype, MPI_Datatype* newtype){

	int blocklength, count, i, countstride;
	MPI_Aint pos, length, max;

	blocklength = k; //Length of '1' blocks

	//count = numprocs/2;	//TODO
	count = badCount(k, numprocs);

	MPI_Aint indices[count];
	int blocklens[count];
	MPI_Datatype oldtypes[count];

	//Get MPI length for sent datatype
	MPI_Type_extent(oldtype, &length);

	//Set max index
	max = numprocs*length;

	//Set start index (see algorithm for details)
	pos = (length)*(rank+k);

	//init data sent count
	countstride=0;

	for(i=0; i<count; i++)
	{
		//Set index of current element
		indices[i] = pos%max;

		//Each element is a block
		blocklens[i] = 1;

		//Oldtype refers to the type sized in length
		oldtypes[i] = oldtype;

		//Count wether we are still on a data to send or not
		++countstride;

		//If we are done for a data block, go to the next block
		if(countstride==blocklength)
		{
			pos+= (length*blocklength) + length; //go to next data block
			countstride=0;	//Reset data block incount
		}

		//else go to the next data inside the block
		else
		{
			pos += length;
		}
	}

	//Create struct
	MPI_Type_create_struct(count, blocklens, indices, oldtypes, newtype);
}

int get_block_to_read_size(int rank, int size, char* file_name, int* buffsize, int * ranks, int * buffs)
{
	int i, k, j;
	char line[255];

	FILE* fh = fopen(file_name, "r");

	i=0;
	j=0;
	k=0;
	*buffsize = 0;
	printf("SIZE %d\n", size);
	while (fgets(line, sizeof(line), fh)) {
		if(i==rank){
			printf("here %d\n", k);
			if(k >= size){
				break;
			}
			else{
				if(line[1] == ':')
					ranks[k] = line[0]-'0';
				else
					ranks[k] = (line[0]-'0')*10 + (line[1]-'0');
				buffs[k] = strlen(line);
				++k;
				*buffsize += strlen(line);
			}
		}
		else{
			k++;
			j+=strlen(line);
			if(k == size)
			{
				k=0;
				i++;
			}
			printf("i %d\n", i);
		}
	}

	fclose(fh);

	i = (*buffsize-1) * rank;
	return j;
}

/**
 * Creates a datatypes for sending via MPI_BOTTOM, 1
 */

void create_read_dt(int rank, int num_proc, int *ranks, int* buffs, char** data, MPI_Datatype* dt, size_t readNum)
 {


	/*
	 * task: Create data structure for reading part
	 */
	assert(data != 0);
	//buffs is the table with the read size

 	//Count variable
 	int i;

 	//Variable for datatype struct almost classic

 	MPI_Aint indices[readNum];
 	int blocklens[readNum];

 	fprintf(stderr, "Rank %d :::::[WRITE] readNum = %zu \n", rank, readNum );
 	MPI_Barrier(MPI_COMM_WORLD);
 	MPI_Datatype oldtypes[readNum];

 	/* Adress originally referencing on data
 	 * data : char** buff in which the read data must be organized by destination rank
 	 */

 	MPI_Aint adress_to_write_in_data_by_element[num_proc];


 	//init adress_to_write_in_data_by_element

 	for(i = 0; i < num_proc; i++){
 		//adress_to_write_in_data_by_element[(rank+i)%num_proc] = (MPI_Aint)data[(rank-i+num_proc)%num_proc];
 		MPI_Get_address(data[(rank-i+num_proc)%num_proc], &adress_to_write_in_data_by_element[(rank+i)%num_proc]);
 	}

 	//double time_count = MPI_Wtime();
 	//Set all classic datatype values
 	for(i = 0; i < readNum; i++){
 		/*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TRICKY PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 		 * Basically what we say here is that the i(th) item that we are going to read goes to the data row corresponding to it's destination rank
 		 *
 		 * indices[i] : The adress for the datatype to which the i(th) item should be written.
 		 * 				Generally this is a relative adress. Like 0, 8, 16.
 		 * 				Here, we use the adress in memory in order to be able to use an array of array without any problem.
 		 * 				Else, the first index could be good
 		 * 				but (data + 8) wouldn't designate &(data[0][8]) if we need to write multiple times in the same row
 		 *
 		 *  adress_to_write_in_data_by_element[i] : This is almost the same as data.
 		 *  		The only difference is that it could be incremented if there are multiple data to write in the same row
 		 *
 		 *  ranks[i] : This designates to what rank the i(th) item should be sent.
 		 *  		Or at least, to what index in data it has to go
 		 *
 		 *  buffs[i] : This is the size of the i(th) item
 		 *
 		 *  blocklens[i] : Nothing special here.
 		 *  		This is the typical blocklens used in MPI struct creation
 		 *
 		 * oldtypes[i] : Nothing special here.
 		 * 			This is the typical oldtypes used in MPI struct creation
 		 */
 		//Set indices
 		//ranks[i] tell the position in adresse to write by elements
 		indices[i] = adress_to_write_in_data_by_element[ranks[i]];



 		//printf("num_proc %d - %d/%d     /    indices: %p   /   buffs: %d    /    ranks: %d\n", size, i, readNum, indices[i], buffs[i], ranks[i]);
 		//Increment position to write for ranks[i]
 		adress_to_write_in_data_by_element[ranks[i]] += buffs[i];

 		/*
 		if (rank == 1 ){
 			fprintf(stderr, "rank %d :::: num_proc %d - %d/%d     /    indices: %p   /   buffs: %d    /    ranks: %d\n",
 							 rank, num_proc, i, readNum, indices[i], buffs[i], ranks[i]);
 		}
 		*/
 		//Set blocklens
 		blocklens[i] = buffs[i];

 		//Set oldtype
 		oldtypes[i] = MPI_CHAR;
 	}


 	for(i = 0; i < readNum; i++){
 		assert (indices[i] != (MPI_Aint)NULL);
 	}


 	//Create struct
 	MPI_Type_create_struct(readNum, blocklens, indices, oldtypes, dt);

 }

int create_send_datatype_for_size(int rank, int size, size_t *num_reads_by_procs, int **dest_size,
		int k, MPI_Datatype* dt, int** recv_index)
{

	/*
	 * task: Create send data type for size
	 */

	int i, j;
	int count = badCount(k, size);
	*recv_index = (int*)malloc(sizeof(int) * count);

	//Variable for datatype struct almost classic
	MPI_Aint indices[count];
	int blocklens[count];
	MPI_Datatype oldtypes[count];

	int *stride = (int *)calloc(k*2, sizeof(int));
	i = 0;
	while (i < (k*2)){
		if (i < k) stride[i] = 0;
		else stride[i] = 1;
		i++;
	}
	int *vect = (int *)calloc(size, sizeof(int));
	i=0;
	while (i < size){
		vect[i] = stride[i%(k*2)];

		i++;
	}

	size_t total=0;
	j = 0;
	for(i = 0; i<size; i++)
	{
		if (vect[i] != 0){
			//MPI_Get_address(dest_size[(i+rank)%size], &indices[j]);
			indices[j] = (MPI_Aint)dest_size[(i+rank)%size];
			(*recv_index)[j] = (i+rank)%size;
			blocklens[j] = (int)(num_reads_by_procs[(i+rank)%size]);
			total += blocklens[j];
			oldtypes[j] = MPI_INT;

			j++;
		}
	}

	//fprintf(stderr, "ran %d ::: indices[0] = %p // dest_offsets[0] = %p \n", rank, indices[0], dest_offsets[(1 + rank)%size]);
	//fprintf(stderr, "ran %d ::: indices[1] = %p // dest_offsets[1] = %p \n", rank, indices[1], dest_offsets[(3 + rank)%size]);
	/*
	for ( i =0; i< count; i++){

		fprintf(stderr, "ran %d ::: blocklens[%d] = %d \n", rank, i, blocklens[i]);
	}
	*/
	free(stride);
	free(vect);

	MPI_Type_create_struct(count, blocklens, indices, oldtypes, dt);

	return count;
}

int create_send_datatype_for_offsets(int rank, int size, size_t *num_reads_by_procs, size_t **dest_offsets,
		int k, MPI_Datatype* dt, int** recv_index)
{

	/*
	 * task: Create send data type for offset
	 */
	int i, j;
	int count = badCount(k, size);
	*recv_index = (int*)malloc(sizeof(int) * count);

	//Variable for datatype struct almost classic
	MPI_Aint indices[count];
	int blocklens[count];
	MPI_Datatype oldtypes[count];

	int *stride = (int *)calloc(k*2, sizeof(int));
	i = 0;
	while (i < (k*2)){
		if (i < k) stride[i] = 0;
		else stride[i] = 1;
		i++;
	}
	int *vect = (int *)calloc(size, sizeof(int));
	i=0;
	while (i < size){
		vect[i] = stride[i%(k*2)];

		i++;
	}

	size_t total=0;
	j = 0;
	for(i = 0; i<size; i++)
	{
		if (vect[i] != 0){
			MPI_Get_address(dest_offsets[(i+rank)%size], &indices[j]);
			//indices[j] = dest_offsets[(i+rank)%size];
			(*recv_index)[j] = (i+rank)%size;
			blocklens[j] = (int)(num_reads_by_procs[(i+rank)%size]);
			total += blocklens[j];
			oldtypes[j] = MPI_LONG_LONG_INT;

			j++;
		}
	}

	//fprintf(stderr, "ran %d ::: indices[0] = %p // dest_offsets[0] = %p \n", rank, indices[0], dest_offsets[(1 + rank)%size]);
	//fprintf(stderr, "ran %d ::: indices[1] = %p // dest_offsets[1] = %p \n", rank, indices[1], dest_offsets[(3 + rank)%size]);
	/*
	for ( i =0; i< count; i++){

		fprintf(stderr, "ran %d ::: blocklens[%d] = %d \n", rank, i, blocklens[i]);
	}
	*/
	free(stride);
	free(vect);

	MPI_Type_create_struct(count, blocklens, indices, oldtypes, dt);

	return count;
}

int create_send_datatype_for_reads(int rank, int size, size_t *buffs, char** data, int k, MPI_Datatype* dt, int** recv_index)
{

	/*
	 * task: Create send data type for reads
	 */


	int i, j;
	int count = badCount(k, size);
	*recv_index = (int*)malloc(sizeof(int) * count);

	//Variable for datatype struct almost classic
	MPI_Aint indices[count];
	int blocklens[count];
	MPI_Datatype oldtypes[count];

	int *stride = (int *)calloc(k*2, sizeof(int));
	i = 0;
	while (i < (k*2)){
		if (i < k) stride[i] = 0;
		else stride[i] = 1;
		i++;
	}
	int *vect = (int *)calloc(size, sizeof(int));
	i=0;
	while (i < size){
		vect[i] = stride[i%(k*2)];

		i++;
	}

	size_t total=0;
	j = 0;
	for(i = 0; i<size; i++)
	{
		if (vect[i] != 0){
			MPI_Get_address(data[(i+rank)%size], &indices[j]);
			(*recv_index)[j] = (i+rank)%size;
			blocklens[j] = (int)buffs[(i+rank)%size];
			total += blocklens[j];
			oldtypes[j] = MPI_CHAR;

			j++;
		}
	}

	free(stride);
	free(vect);

	MPI_Type_create_struct(count, blocklens, indices, oldtypes, dt);

	return count;
}

size_t get_send_size(int rank, int size, size_t* buffs, size_t** send_size, int count, int k)
{
	int i, j;

	int *stride = (int *)calloc(k*2, sizeof(int));
	i = 0;
	while (i < (k*2)){
		if (i < k) stride[i] = 0;
		else stride[i] = 1;
		i++;
	}
	int *vect = (int *)calloc(size, sizeof(int));
	i=0;
	while (i < size){
		vect[i] = stride[i%(k*2)];
		i++;
	}

	size_t total=0;
	j=0;
	for(i=0; i<size; i++)
	{
		if(vect[i] != 0)
		{
			total += buffs[(i+rank)%size];
			(*send_size)[j] = buffs[(i+rank)%size];

			j++;
		}
	}

	free(stride);
	free(vect);

	return total;
}

int badCount(int k, int n)
{
	int c = ((int)(n/(k*2) ))*k;
	//printf("\tc : %d\n", c);
	if(k>1)
		if((int)((n%(k*2) )-k) > 0)
			c += (int)((n%(k*2) )-k);
	//printf("\tc (encore) : %d\n", c);
	return c;
}
