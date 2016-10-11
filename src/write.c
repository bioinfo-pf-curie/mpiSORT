#include "write.h"
#include "bgzf.c"
#include "bgzf.h"
#include "mpiSort_utils.h"

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

void read_data_for_writing(int rank, int num_proc, size_t local_readNum, char *file_name,
		size_t *number_of_reads_by_procs, size_t *buffs_by_procs, char *** data,
		int *new_rank, int *new_size, size_t *new_offset, MPI_File in, MPI_Info finfo, MPI_Comm COMM_WORLD)
{

	/*
	 * task: IN read_data_for_writing
	 */

	size_t k;
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

	//t = MPI_Wtime();
	MPI_File_set_view(in, 0, MPI_CHAR, dt_view, "native", finfo);
	MPI_File_read(in, MPI_BOTTOM, 1, dt_data, MPI_STATUS_IGNORE);
	//fprintf(stderr, "%d ::::: [read_data_for_writing] Time in MPI_File_read %f sec\n", rank, MPI_Wtime() - t);

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

void bruck_reads(int rank, int num_proc, size_t * buffs_by_procs, char** data2)
{
	MPI_Comm comm = COMM_WORLD;

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
				data_offsets[recv_index[m]] = (size_t *)malloc(sizeof(size_t)*(recv_size_by_proc[m]));
				// TODO Problem when writing
				// valgrind reports invalid write
				//data_offsets[recv_index[m]][0] = 0;
			}
				//data_offsets[recv_index[m]] = realloc(data_offsets[recv_index[m]], sizeof(size_t)*(recv_size_by_proc[m]));

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
	MPI_Comm comm = COMM_WORLD;

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
				data_size[recv_index[m]] = NULL;
				data_size[recv_index[m]] = (int *)malloc(sizeof(int)*(recv_size_by_proc[m]));
				// TODO Problem when writing
				// valgrind reports invalid write
				//data_size[recv_index[m]][0] = 0;

				//data_size[recv_index[m]] = realloc(data_size[recv_index[m]], sizeof(int)*(recv_size_by_proc[m]) );

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
		char *data,
		size_t start_offset_in_file){


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

	int *all_rank_to_send=NULL;
	size_t *all_offset_dest_file_to_send=NULL;
	int *all_read_size_to_send=NULL;
	size_t *all_offset_source_file_to_send=NULL;

	/*
	 * phase 1 variables
	 */

	int *all_read_size_phase1_to_send=NULL;
	int *all_rank_phase1_to_send=NULL;
	size_t* all_offset_source_file_phase1_to_send=NULL;
	size_t *new_offset_dest_phase1 = (size_t *) malloc(local_readNum* sizeof(size_t));
	size_t *new_offset_source_phase1 = (size_t *) malloc(local_readNum* sizeof(size_t));
	size_t *new_offset_source_sorted_phase1 = (size_t *) malloc(local_readNum* sizeof(size_t));
	int *new_read_size_phase1 = (int *) malloc(local_readNum* sizeof(int));
	int *new_rank_phase1 = (int *) malloc(local_readNum* sizeof(int));
	size_t *new_offset_dest_index_phase1 = (size_t *) malloc(local_readNum* sizeof(size_t));
	size_t *pbs_local_dest_offset=NULL;
	size_t *pbs_local_dest_offset_index=NULL;
	size_t *all_offset_dest_sorted_index_phase1=NULL;
	size_t *all_offset_dest_file_to_send_phase1=NULL;

	/*
	 * phase 2 variables
	 */

	size_t *all_offset_source_file_to_send_phase2=NULL;
	int *all_read_size_phase2_to_send=NULL;
	size_t *all_offset_dest_file_phase2_to_send=NULL;
	int *all_rank_phase2_to_send=NULL;
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
		fprintf(stderr, "Rank %d :::::[WRITE][PHASE 1] SORT OFFSET DESTINATION \n", rank);

	size_t *num_reads_per_jobs = (size_t *) malloc(total_num_proc * sizeof(size_t));

	MPI_Gather(&local_readNum, 1, MPI_LONG_LONG_INT, &num_reads_per_jobs[rank - master_job_phase_2], 1,
			MPI_LONG_LONG_INT, master_job_phase_1 , COMM_WORLD);


		//vector of all offset in destination / source file / size of reads / ranks
		size_t *all_offset_dest_file_phase1=NULL;
		size_t *all_offset_source_file_phase1=NULL;
		int *all_read_size_phase1=NULL;
		int *all_rank_phase1=NULL;

		size_t total_num_read_phase1 = 0;
		MPI_Reduce(&local_readNum, &total_num_read_phase1, 1, MPI_LONG_LONG_INT, MPI_SUM, master_job_phase_1, COMM_WORLD);

		if (rank == master_job_phase_1)
			assert(total_num_read_phase1 == total_num_read);

		if (rank == master_job_phase_1){
			all_offset_dest_file_phase1 = (size_t *) malloc (total_num_read_phase1 * sizeof(size_t));
			all_offset_source_file_phase1 = (size_t *) malloc (total_num_read_phase1 * sizeof(size_t));
			all_read_size_phase1 = (int *) malloc (total_num_read_phase1 * sizeof(int));
			all_rank_phase1 = (int *) malloc (total_num_read_phase1 * sizeof(int));

			size_t k =0;

			for (k = 0; k < total_num_read_phase1; k++){
				all_read_size_phase1[k] = 0;
				all_offset_dest_file_phase1[k] = 0;
				all_offset_source_file_phase1[k] = 0;
				all_rank_phase1[k] = 0;
			}
		}

		/*
		 * Phase 1: master_1 defines the following two vectors
		 */
		// vector of number of read per jobs
		size_t *num_reads_per_jobs_phase1 = (size_t *) malloc(total_num_proc* sizeof(size_t));
		// vector of index that contains the cumulative sum of the number of reads
		size_t *start_num_reads_per_jobs_phase1 = (size_t *) malloc((total_num_proc + 1)*sizeof(size_t));

		/*
		 * Phase 2: master_2 receives all local_readNum and adds it to a local vector
		 */
		MPI_Gather(&local_readNum, 1, MPI_LONG_LONG_INT, &num_reads_per_jobs_phase1[rank - master_job_phase_1],
				1, MPI_LONG_LONG_INT, master_job_phase_1, COMM_WORLD);

		if (rank == master_job_phase_1){

			start_num_reads_per_jobs_phase1[0] = 0;

			for (k = 1; k < (total_num_proc +1); k++){
				start_num_reads_per_jobs_phase1[k] = num_reads_per_jobs_phase1[k-1];
			}

			for (k = 1; k < total_num_proc; k++){
				size_t tmp = start_num_reads_per_jobs_phase1[k - 1];
				size_t tmp2 = start_num_reads_per_jobs_phase1[k];
				start_num_reads_per_jobs_phase1[k] = tmp + tmp2;
			}
		}

		if (rank == master_job_phase_1){

			MPI_Status status;
			//we copy the first elements in
			int k=0;
			size_t st = start_num_reads_per_jobs_phase1[master_job_phase_1];

			for (k = 0; k < num_reads_per_jobs_phase1[master_job_phase_1]; k++){

				all_offset_dest_file_phase1[st] = new_offset_dest[k];
				all_offset_source_file_phase1[st] = new_offset_source[k];
				all_read_size_phase1[st] = new_read_size[k];
				all_rank_phase1[st] = new_rank[k];
				st++;
			}

			for(j = 0; j < total_num_proc; j++){

				if(j != master_job_phase_2){

					// first we care for ranks
					int *temp_buf =(int *) malloc(num_reads_per_jobs_phase1[j]* sizeof(int));
					int *temp_buf1 =(int *) malloc(num_reads_per_jobs_phase1[j]* sizeof(int));
					size_t *temp_buf2 =(size_t *) malloc(num_reads_per_jobs_phase1[j]* sizeof(size_t));
					size_t *temp_buf3 =(size_t *) malloc(num_reads_per_jobs_phase1[j]* sizeof(size_t));

					MPI_Recv(temp_buf, num_reads_per_jobs_phase1[j], MPI_INT, j, 0, COMM_WORLD, &status);
					MPI_Recv(temp_buf1, num_reads_per_jobs_phase1[j], MPI_INT, j, 1, COMM_WORLD, &status);
					MPI_Recv(temp_buf2, num_reads_per_jobs_phase1[j], MPI_LONG_LONG_INT, j, 2, COMM_WORLD , &status);
					MPI_Recv(temp_buf3, num_reads_per_jobs_phase1[j], MPI_LONG_LONG_INT, j, 3, COMM_WORLD, &status);

					st=0;
					size_t st = start_num_reads_per_jobs_phase1[j];

					for (k = 0; k < num_reads_per_jobs_phase1[j]; k++){

						all_rank_phase1[st] = temp_buf[k];
						all_read_size_phase1[st] = temp_buf1[k];
						all_offset_source_file_phase1[st] = temp_buf2[k];
						all_offset_dest_file_phase1[st] = temp_buf3[k];
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

		int dimensions = 1;
		while (dimensions <= total_num_proc)
			dimensions <<= 1;

		dimensions >>= 1;

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

		while((pbs_num_offsets_to_recieve_left_phase1 * dimensions) > pbs_num_offsets_to_recieve_phase1){
			dimensions >>= 1;
			pbs_num_offsets_to_recieve_phase1 = total_num_read_phase1/dimensions;
			pbs_num_offsets_to_recieve_left_phase1  = total_num_read_phase1 - pbs_num_offsets_to_recieve_phase1*dimensions;
		}


		// we compute a vector of size dimensions which contain the number
		// of reads to send
		size_t *pbs_local_num_read_per_job_phase1 = (size_t *)malloc(dimensions * sizeof(size_t));

		// each job create a vector with with the length
		// of the offset vector for each job in [0, dimension[
		for (j = 0; j < dimensions ; j++){
			// we add pbs_num_offsets_to_recieve_left
			// because we want the vector to have the same size
			pbs_local_num_read_per_job_phase1[j] = pbs_num_offsets_to_recieve_phase1 + pbs_num_offsets_to_recieve_left_phase1;
		}
		size_t *pbs_start_num_offset_per_jobs_phase1 = (size_t *) malloc((dimensions + 1) * sizeof(size_t));

		// the master job compute the start index of the element
		// to dispatch
		if (rank == master_job_phase_1){

			pbs_start_num_offset_per_jobs_phase1[0] = 0;

			for (k = 1; k < (dimensions +1); k++){
				pbs_start_num_offset_per_jobs_phase1[k] = pbs_local_num_read_per_job_phase1[k-1];
			}
			for (k = 1; k < dimensions; k++){
				size_t tmp = pbs_start_num_offset_per_jobs_phase1[k - 1];
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
			pbs_local_dest_offset = (size_t *)malloc(sizeof(size_t) * pbs_local_num_read_per_job_phase1[rank]);
			//now the master send

			if ( rank != master_job_phase_1 ){
					MPI_Recv(pbs_local_dest_offset, pbs_local_num_read_per_job_phase1[rank], MPI_LONG_LONG_INT,
							master_job_phase_1, 0, COMM_WORLD, &status);

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
						MPI_Send(&all_offset_dest_file_phase1[pbs_start_num_offset_per_jobs_phase1[j]],
							pbs_local_num_read_per_job_phase1[j], MPI_LONG_LONG_INT, j, 0, COMM_WORLD);
					}
				}
			}

			// we build pbs_local_dest_offset_index
			pbs_local_dest_offset_index = (size_t *)malloc(pbs_local_num_read_per_job_phase1[rank]*sizeof(size_t));

			for (j = 0; j < pbs_local_num_read_per_job_phase1[rank]; j++){

				if (rank == master_job_phase_1){
					pbs_local_dest_offset_index[j] = j + pbs_local_num_read_per_job_phase1[rank]*rank;
				}
				else{
					pbs_local_dest_offset_index[j] = j + pbs_local_num_read_per_job_phase1[rank]*rank - (rank*pbs_num_offsets_to_recieve_left_phase1);
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
				all_read_size_phase1_to_send = (int*)malloc(total_num_read_phase1*sizeof(int));
				all_offset_source_file_phase1_to_send = (size_t*)malloc(total_num_read_phase1*sizeof(size_t));
				all_rank_phase1_to_send = (int*)malloc(total_num_read_phase1*sizeof(int));

				for(j = 0; j < total_num_read_phase1; j++){

					all_offset_source_file_phase1_to_send[j] = all_offset_source_file_phase1[all_offset_dest_sorted_index_phase1[j]];
					all_read_size_phase1_to_send[j] = all_read_size_phase1[all_offset_dest_sorted_index_phase1[j]];

				}


				//now we change the rank
				// we initialize all_offset_rank_to_send constains
				// the rank of the sorted read
				/*
				 * improvement
				 *
				 *
				 */
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
			free(pbs_local_num_read_per_job_phase1);
			free(pbs_start_num_offset_per_jobs_phase1);
		}

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

				new_offset_dest_phase1[k] = all_offset_dest_file_to_send_phase1[ind];
				new_offset_source_phase1[k] = all_offset_source_file_phase1_to_send[ind];
				new_read_size_phase1[k] = all_read_size_phase1_to_send[ind];
				new_rank_phase1[k] = all_rank_phase1_to_send[ind];
				ind++;
			}

			for(j = 0; j < total_num_proc; j++){

				if (j != master_job_phase_1){

					MPI_Send(&all_offset_dest_file_to_send_phase1[start_num_reads_per_jobs_phase1[j]],
							num_reads_per_jobs_phase1[j], MPI_LONG_LONG_INT, j, 0, COMM_WORLD);

					MPI_Send(&all_offset_source_file_phase1_to_send[start_num_reads_per_jobs_phase1[j]],
							num_reads_per_jobs_phase1[j], MPI_LONG_LONG_INT, j, 1, COMM_WORLD);

					MPI_Send(&all_read_size_phase1_to_send[start_num_reads_per_jobs_phase1[j]],
							num_reads_per_jobs_phase1[j], MPI_INT, j, 2, COMM_WORLD);

					MPI_Send(&all_rank_phase1_to_send[start_num_reads_per_jobs_phase1[j]],
							num_reads_per_jobs_phase1[j], MPI_INT, j, 3, COMM_WORLD);

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
	MPI_Gather(&local_readNum, 1, MPI_LONG_LONG_INT, &num_reads_per_jobs[rank - master_job_phase_2], 1,
			MPI_LONG_LONG_INT, master_job_phase_2 , COMM_WORLD);
	if (rank == master_job_phase_2)
		fprintf(stderr, "Rank %d :::::[WRITE][PHASE 2] SORT OFFSET SOURCES \n", rank);



	//vector of all offset in destination / source file / size of reads / ranks
	size_t *all_offset_dest_file_phase2=NULL;
	size_t *all_offset_source_file_phase2=NULL;
	int *all_read_size_phase2=NULL;
	int *all_rank_phase2=NULL;

	/*
	 * Phase 2: master_2 gets the total of reads
	 * Improvement: see if the total num read is the same of the phase1
	 */


	size_t total_num_read_phase2 = 0;
	MPI_Reduce(&local_readNum, &total_num_read_phase2, 1, MPI_LONG_LONG_INT, MPI_SUM, master_job_phase_2, COMM_WORLD);

	if (rank == master_job_phase_2)
		assert(total_num_read_phase2 == total_num_read);

	if (rank == master_job_phase_2){
		all_offset_dest_file_phase2 = (size_t *) malloc (total_num_read_phase2 * sizeof(size_t));
		all_offset_source_file_phase2 = (size_t *) malloc (total_num_read_phase2 * sizeof(size_t));
		all_read_size_phase2 = (int *) malloc (total_num_read_phase2 * sizeof(int));
		all_rank_phase2 = (int *) malloc (total_num_read_phase2 * sizeof(int));

		size_t k =0;

		for (k = 0; k < total_num_read_phase2; k++){
			all_read_size_phase2[k] = 0;
			all_offset_dest_file_phase2[k] = 0;
			all_offset_source_file_phase2[k] = 0;
			all_rank_phase2[k] = 0;
		}
	}

	// vector of number of read per jobs
	size_t *num_reads_per_jobs_phase2 = (size_t *) malloc(total_num_proc* sizeof(size_t));
	// vector of index that contains the cumulative sum of the number of reads
	size_t *start_num_reads_per_jobs_phase2 = (size_t *) malloc((total_num_proc + 1)*sizeof(size_t));

	/*
	 * Phase 2: master_2 receives all local_readNum and adds it to a local vector
	 */
	MPI_Gather(&local_readNum, 1, MPI_LONG_LONG_INT, &num_reads_per_jobs_phase2[rank - master_job_phase_2],
			1, MPI_LONG_LONG_INT, master_job_phase_2, COMM_WORLD);


	if (rank == master_job_phase_2){

		start_num_reads_per_jobs_phase2[0] = 0;

		for (k = 1; k < (total_num_proc +1); k++){
			start_num_reads_per_jobs_phase2[k] = num_reads_per_jobs_phase2[k-1];
		}

		for (k = 1; k < total_num_proc; k++){
			size_t tmp = start_num_reads_per_jobs_phase2[k - 1];
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

			all_offset_dest_file_phase2[st] = new_offset_dest_phase1[k];
			all_offset_source_file_phase2[st] = new_offset_source_phase1[k];
			all_read_size_phase2[st] = new_read_size_phase1[k];
			all_rank_phase2[st] = new_rank_phase1[k];
			st++;
		}

		for(j = 0; j < total_num_proc; j++){

			if(j != master_job_phase_2){

				// first we care for ranks
				int *temp_buf =(int *) malloc(num_reads_per_jobs_phase2[j]* sizeof(int));
				int *temp_buf1 =(int *) malloc(num_reads_per_jobs_phase2[j]* sizeof(int));
				size_t *temp_buf2 =(size_t *) malloc(num_reads_per_jobs_phase2[j]* sizeof(size_t));
				size_t *temp_buf3 =(size_t *) malloc(num_reads_per_jobs_phase2[j]* sizeof(size_t));

				MPI_Recv(temp_buf, num_reads_per_jobs_phase2[j], MPI_INT, j, 0, COMM_WORLD, &status);
				MPI_Recv(temp_buf1, num_reads_per_jobs_phase2[j], MPI_INT, j, 1, COMM_WORLD, &status);
				MPI_Recv(temp_buf2, num_reads_per_jobs_phase2[j], MPI_LONG_LONG_INT, j, 2, COMM_WORLD , &status);
				MPI_Recv(temp_buf3, num_reads_per_jobs_phase2[j], MPI_LONG_LONG_INT, j, 3, COMM_WORLD, &status);

				st=0;
				size_t st = start_num_reads_per_jobs_phase2[j];

				for (k = 0; k < num_reads_per_jobs_phase2[j]; k++){

					all_rank_phase2[st] = temp_buf[k];
					all_read_size_phase2[st] = temp_buf1[k];
					all_offset_source_file_phase2[st] = temp_buf2[k];
					all_offset_dest_file_phase2[st] = temp_buf3[k];
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
			assert ( all_offset_dest_file_phase2[j] != 0 );
			assert ( all_offset_source_file_phase2[j] != 0 );
			assert ( all_read_size_phase2[j] != 0 );
			assert ( all_rank_phase2[j] <= total_num_proc);
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

	dimensions = 1;
	while (dimensions <= total_num_proc)
		dimensions <<= 1;

	dimensions >>= 1;

	// the master rank compute the number of
	// dimension is the number of processors where we
	// perform the bitonic sort
	// int dimensions = (int)(log2(num_processes));
	// find next ( must be greater) power, and go one back

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
	size_t *pbs_local_num_read_per_job_phase2 = (size_t *)malloc(dimensions * sizeof(size_t));

	// each job create a vector with with the length
	// of the offset vector for each job in [0, dimension[
	for (j = 0; j < dimensions ; j++){
		// we add pbs_num_offsets_to_recieve_left
		// because we want the vector to have the same size
		pbs_local_num_read_per_job_phase2[j] = pbs_num_offsets_to_recieve + pbs_num_offsets_to_recieve_left;
	}


	size_t *pbs_start_num_offset_per_jobs_phase2 = (size_t *) malloc((dimensions + 1) * sizeof(size_t));

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
		pbs_local_source_offset = (size_t *)malloc(sizeof(size_t) * pbs_local_num_read_per_job_phase2[rank]);
		//now the master send

		if ( rank != master_job_phase_2 ){
				MPI_Recv(pbs_local_source_offset, pbs_local_num_read_per_job_phase2[rank], MPI_LONG_LONG_INT,
						master_job_phase_2, 0, COMM_WORLD, &status);

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
						pbs_local_num_read_per_job_phase2[j], MPI_LONG_LONG_INT, j, 0, COMM_WORLD);
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

		for(j = 0; j < pbs_local_num_read_per_job_phase2[rank]; j++){
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
						start_index);

		if (rank == master_job_phase_2){
			fprintf(stderr, "Rank %d :::::[WRITE][PHASE 3] Time for gathering source offset and index %f\n", rank, MPI_Wtime() - time_count);


			// now we apply the new index to all the
			all_read_size_phase2_to_send = (int*)malloc(total_num_read_phase2*sizeof(int));
			all_offset_dest_file_phase2_to_send = (size_t*)malloc(total_num_read_phase2*sizeof(size_t));
			all_rank_phase2_to_send = (int*)malloc(total_num_read_phase2*sizeof(int));

			for(j = 0; j < total_num_read_phase2; j++){

				all_offset_dest_file_phase2_to_send[j] = all_offset_dest_file_phase2[all_offset_source_sorted_index_phase2[j]];
				all_read_size_phase2_to_send[j] = all_read_size_phase2[all_offset_source_sorted_index_phase2[j]];
				all_rank_phase2_to_send[j] = all_rank_phase2[all_offset_source_sorted_index_phase2[j]];
			}

			if (rank == master_job_phase_2){
				for(j = 0; j < total_num_read_phase2; j++){
					assert ( all_offset_dest_file_phase2_to_send[j] != 0 );
					assert ( all_read_size_phase2_to_send[j] != 0 );
					assert ( all_rank_phase2_to_send[j] <= total_num_proc );
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
	size_t *buffs_by_procs = (size_t*)calloc( total_num_proc, sizeof(size_t));

	time_count = MPI_Wtime();
	int num_proc = total_num_proc;


	/*
	 *
	read_data_for_writing(
				rank,
				num_proc,
				local_readNum,
				file_name,
				number_of_reads_by_procs,
				buffs_by_procs,
				&data2,
				new_rank_phase2,
				new_read_size_phase2,
				new_offset_source_sorted_phase2,
				in,
				finfo,
				COMM_WORLD);
	 */

	time_count = MPI_Wtime();


	/*
	 * first we pack data in data_pack1
	 * to have consecutive reads of the
	 * chromosoms
	 */

	int m;
	size_t i;
	size_t new_data_sz = 0;

	MPI_Datatype dt_data0;
	char *data_pack;

	//we compute the size of data_pack
	for (k = 0; k < local_readNum; k++){
		new_data_sz += new_read_size_phase2[k];
	}

	data_pack = malloc(new_data_sz +1);
	data_pack[new_data_sz] = 0;

	//we create a new offset vector
	//wich is the offset in the data
	size_t *offset_in_data = (size_t *) malloc(local_readNum* sizeof(size_t));

	//we compute the new offset of reads in data buffer
	//we remove the start offset in the file

	for (k=0; k < local_readNum;k++){
		offset_in_data[k] = new_offset_source_sorted_phase2[k] - start_offset_in_file;
	}

	MPI_Type_hindexed(local_readNum, &new_read_size_phase2[0], (MPI_Aint*)offset_in_data, MPI_CHAR, &dt_data0);
	MPI_Type_commit(&dt_data0);

	int position=0;

	MPI_Pack(data, 1, dt_data0, data_pack, new_data_sz, &position, COMM_WORLD);
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

	//Variable for datatype struct almost classic
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
	 int pos=0;
	 MPI_Unpack(data_pack, new_data_sz, &pos, MPI_BOTTOM, 1, dt_data, COMM_WORLD);

	/*
	 if (rank == master_job_phase_2){
		int j =0;
		for (j = 0; j < num_proc; j++){
			fprintf(stderr, "rank %d ::::: [read_data_for_writing] data2 = \n", rank);
			for (m = 0; m < buffs_by_procs[j]; m++){
				fprintf(stderr, "%c", data2[j][m]);
			}
			fprintf(stderr, "\n");
		}
	 }
	*/
	MPI_Type_free(&dt_data);
	MPI_Type_free(&dt_data0);
	free(data_pack);
	free(blocklens);
	free(indices);
	free(oldtypes);
	free(new_offset_source_sorted_phase2);
	if (rank == master_job_phase_2)
		fprintf(stderr, "Rank %d :::::[WRITE][PACK] Time in packing data for bruck %f seconds\n", rank, MPI_Wtime() - time_count);

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

	if (rank == master_job_phase_2)
		fprintf(stderr, "Rank %d :::::[WRITE][BRUCK] Time for Bruck phases %f seconds\n", rank, MPI_Wtime() - time_count);

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

	//int i=0;
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
	char_buff_uncompressed[y_message_sz[predecessor]] = 0;

	//now we send data
	MPI_Sendrecv(MPI_BOTTOM, 1, Datatype_Read_to_write, successor, 0, char_buff_uncompressed, y_message_sz[predecessor],
            MPI_CHAR, predecessor, 0, COMM_WORLD,  &status);

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

	free(char_buff_uncompressed);

	MPI_Barrier(COMM_WORLD);

	//we trade the blocks
	//in this phase the predeccor become the succesor
	size_t compressed_sz_to_send = compressed_size;
	size_t compressed_sz_to_recv = 0;

	// copy x into the correct location in y
	int predecessor_back = ( rank + 1 ) % num_proc;
	int successor_back = ( rank - 1 + num_proc ) % num_proc;

	// first we size the size of the buffer to send
	MPI_Send( &compressed_sz_to_send, blocksize , MPI_LONG_LONG_INT, successor_back, 0, COMM_WORLD);
	MPI_Recv( &compressed_sz_to_recv, blocksize , MPI_LONG_LONG_INT, predecessor_back, 0, COMM_WORLD, &status);

	//we create a buffer for recieved data
	uint8_t *buff_compressed = malloc(compressed_sz_to_recv * sizeof(uint8_t));
	//char *buff_compressed = malloc(compressed_sz_to_recv * sizeof(char));
	//now we send data
	MPI_Sendrecv(compressed_buff, compressed_sz_to_send, MPI_UNSIGNED_CHAR, successor_back, 0,
			buff_compressed, compressed_sz_to_recv, MPI_UNSIGNED_CHAR, predecessor_back, 0, COMM_WORLD,  &status);

	free(compressed_buff);

	MPI_Barrier(COMM_WORLD);
	size_t compSize = compressed_sz_to_recv;


	/*
	 * Now we write results of compression
	 */
	MPI_Barrier(COMM_WORLD);
	size_t write_offset = 0;

	MPI_Offset * y = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
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


	/*
	 * now we compute offset where to write
	 */

	// we create the path where to write for collective write
	path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
	sprintf(path, "%s/%s.bam", output_dir, chrName);

	if(!rank)
		fprintf(stderr, "rank %d :::: Opening the file %s \n", rank, path );


	//task FINE TUNING FINFO FOR WRITING OPERATIONS
	MPI_Info_set(finfo,"striping_factor","128");
	MPI_Info_set(finfo,"striping_unit","1610612736"); //1G striping

	MPI_Info_set(finfo,"nb_proc","128");
	MPI_Info_set(finfo,"cb_nodes","128");
	MPI_Info_set(finfo,"cb_block_size","1610612736"); /* 4194304 MBytes - should match FS block size */
	MPI_Info_set(finfo,"cb_buffer_size","1610612736"); /* 128 MBytes (Optional) */


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

	if (rank == master_job_phase_2 ) {
		fprintf(stderr, "Proc rank %d ::: we write the header \n", rank);
		MPI_File_write(out, compressed_header, compressed_size_header, MPI_BYTE, MPI_STATUS_IGNORE);
	}
	free(compressed_header);

	MPI_File_set_view(out, write_offset, MPI_BYTE, MPI_BYTE, "native", finfo);
	MPI_File_write_all(out, buff_compressed, (size_t)compSize, MPI_BYTE, &status);



	//task FINE TUNING FINFO BACK TO READING OPERATIONS
	MPI_Info_set(finfo,"striping_factor","128");
	MPI_Info_set(finfo,"striping_unit","2684354560"); //1G striping

	MPI_Info_set(finfo,"nb_proc","128");
	MPI_Info_set(finfo,"cb_nodes","128");
	MPI_Info_set(finfo,"cb_block_size","2684354560"); /* 4194304 MBytes - should match FS block size */
	MPI_Info_set(finfo,"cb_buffer_size","2684354560"); /* 128 MBytes (Optional) */

	free(buff_compressed);

	if (rank == master_job_phase_2)
		fprintf(stderr, "Rank %d :::::[WRITE] Time for chromosome %s writing %f seconds\n",
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
	free(data2);
	free(data_table);
	free(data_reads_to_sort);
	free(reads_address_sorted);
	free(offsets_sorted);
	free(new_read_size_sorted_phase3);
	free(data_size_to_sort);

}

void writeSam_discordant_and_unmapped(int split_rank, char* output_dir, char* header, size_t local_readNum, char* chrName, Read* chr,
		int split_size, MPI_Comm split_comm, char *file_name, MPI_File in, MPI_Info finfo, int compression_level,
		char* data, size_t start_offset_in_file){

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

	int master_job = 0;
	double start_phase2, finish_phase2;

	//we initialize offset source and size_source
	for(j = 0; j < local_readNum; j++){
		 size_source[j] = 0;
		 offset_source[j] = 0;
	}

	int *all_read_size=NULL;
	size_t *all_offset_source_file=NULL;
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
	MPI_Reduce(&local_readNum, &total_num_read, 1, MPI_LONG_LONG_INT, MPI_SUM, master_job, split_comm);

	MPI_Bcast(&total_num_read, 1, MPI_LONG_LONG_INT, master_job, split_comm );

	if (split_rank == master_job){
		all_offset_source_file = (size_t *) malloc (total_num_read * sizeof(size_t));
		all_read_size = (int *) malloc (total_num_read * sizeof(int));
	}


	// vector of number of read per jobs
	size_t *num_reads_per_jobs = (size_t *) malloc(split_size* sizeof(size_t));
	// vector of index that contains the cumulative sum of the number of reads
	size_t *start_num_reads_per_jobs = (size_t *) malloc((split_size + 1)*sizeof(size_t));

	// job 1 recieves the number
	// of reads of each rank
	// and put it  in a vector

	MPI_Gather(&local_readNum, 1, MPI_LONG_LONG_INT, &num_reads_per_jobs[split_rank - master_job], 1,
			MPI_LONG_LONG_INT, master_job , split_comm);

	if (split_rank == master_job){

		start_num_reads_per_jobs[0] = 0;

		for (k = 1; k < (split_size +1); k++){
			start_num_reads_per_jobs[k] = num_reads_per_jobs[k-1];
		}
		for (k = 1; k < split_size; k++){
			size_t tmp = start_num_reads_per_jobs[k - 1];
			size_t tmp2 = start_num_reads_per_jobs[k];
			start_num_reads_per_jobs[k] = tmp + tmp2;
		}
	}

	if (split_rank == master_job){

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
		 *  	num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, MPI_split_comm, &status);
		 *
		 */

		//fprintf(stderr, "Rank %d ::::: Phase 2 :::: recieve from other job", rank);
		for(j = 0; j < split_size; j++){

			if ( j != master_job ){

				size_t *temp_buf =(size_t *) malloc(num_reads_per_jobs[j]* sizeof(size_t));
				MPI_Recv(temp_buf, num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, split_comm, &status);

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
		MPI_Send(new_offset_source, local_readNum, MPI_LONG_LONG_INT, master_job,  0, split_comm);
	}


	//we create vector with all the reads sizes
	//attached to the offset in dest file

	if (split_rank == master_job){

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
		 *  	num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, MPI_split_comm, &status);
		 *
		 */

		for(j = 0; j < split_size; j++){

			if (j != master_job){

				int *temp_buf =(int *) malloc(num_reads_per_jobs[j]* sizeof(int));
				MPI_Recv(temp_buf, num_reads_per_jobs[j], MPI_INT, j, 0, split_comm, &status);

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
		MPI_Send(new_read_size, local_readNum, MPI_INT, master_job,  0, split_comm);
	}

	/**************************/
	// We free some variable
	/**************************/
	free(new_read_size);
	free(new_offset_dest);
	free(new_offset_source);
	free(new_rank);


	size_t *all_reads_offset_source_index = NULL;

	//task: Phase two: Sorting all offset sources with bitonic
	int chosen_split_rank= 0;
	int dimensions = 1;
	while (dimensions <= split_size)
		dimensions <<= 1;

	dimensions >>= 1;

	// we compute the length of the vector to recieve
	// the we split all_offset_dest_index_phase1 among the dimension processors
	size_t pbs_num_coordinates_to_recieve = total_num_read/dimensions;
	// for the last processor will recieve the rest of the division
	size_t pbs_num_coordinates_to_recieve_left  =  total_num_read - pbs_num_coordinates_to_recieve * dimensions;

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
	//pbs_local_num_read_per_job[dimensions - 1] += pbs_num_coordinates_to_recieve_left;
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

	size_t *pbs_local_offset_source_file = NULL;
	MPI_Barrier(split_comm);

	if (split_rank < dimensions){

		// pbs_local_reads_coordinates is a table containing the unsorted
		// reference coordinates
		pbs_local_offset_source_file = (size_t *)malloc(sizeof(size_t) * pbs_local_num_read_per_job[split_rank]);
		//pbs_local_reads_sizes = (int *)malloc(sizeof(int) * pbs_local_num_read_per_job[split_rank]);

		//now the master send
		//fprintf(stderr, "Rank %d :::::[WRITE] dispatch of reads coordinates for bitonic dimensions = %d \n", split_rank, dimensions);
		time_count = MPI_Wtime();
		if ( split_rank != chosen_split_rank ){

			MPI_Status status;
			//fprintf(stderr, "Rank %d :::::[WRITE] MPI_Recv from %d\n", rank, master_job_phase_1);
			MPI_Recv(pbs_local_offset_source_file, pbs_local_num_read_per_job[split_rank],
					MPI_LONG_LONG_INT, chosen_split_rank, 0, split_comm, &status);

		}
		else {
			//first we copy the data from the master job
			size_t ind = pbs_start_num_coordinates_per_jobs[chosen_split_rank];
			int k = 0;

			for (k = 0; k < pbs_local_num_read_per_job[chosen_split_rank]; k++){
				pbs_local_offset_source_file[k] = all_offset_source_file[ind];
				ind++;
			}

			for(j = 0; j < dimensions; j++){
				if (j != chosen_split_rank){
					MPI_Send(&all_offset_source_file[pbs_start_num_coordinates_per_jobs[j]],
							pbs_local_num_read_per_job[j], MPI_LONG_LONG_INT, j, chosen_split_rank, split_comm);
				}
			}
		}
		if (split_rank == chosen_split_rank)
			fprintf(stderr,	"rank %d :::::[mpiSort] Time to dispatch all_read_size  %f seconds\n", split_rank, MPI_Wtime() - time_count);




		// we build pbs_local_reads_coordinates_index
		size_t *pbs_global_reads_offset_source_index = (size_t *)malloc(pbs_local_num_read_per_job[split_rank]*sizeof(size_t));

		for (j = 0; j < pbs_local_num_read_per_job[split_rank]; j++){

			if (split_rank == chosen_split_rank){
				pbs_global_reads_offset_source_index[j] = j + pbs_local_num_read_per_job[split_rank] * split_rank;
			}
			else{
				pbs_global_reads_offset_source_index[j] = j + pbs_local_num_read_per_job[split_rank] * split_rank -
						(split_rank * pbs_num_coordinates_to_recieve_left);
			}
		}

		if (split_rank == chosen_split_rank)
					fprintf(stderr,	"rank %d :::::[mpiSort] We call bitonic with dimensions = %d \n", split_rank, dimensions);

		// now each rank from [0, dimension[
		// is going to bitonic sort

		time_count = MPI_Wtime();
		ParallelBitonicSort(split_comm, split_rank, dimensions, pbs_local_offset_source_file,
				pbs_global_reads_offset_source_index, pbs_local_num_read_per_job[split_rank],
				pbs_num_coordinates_to_recieve_left);

		if (split_rank == chosen_split_rank)
			fprintf(stderr,	"rank %d :::::[mpiSort] Time in parallel bitonic sort  %f seconds\n", split_rank, MPI_Wtime() - time_count);

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
		while (pbs_local_offset_source_file[start_index] == 0){
			start_index++;
		}

		free(pbs_local_offset_source_file); //ok

		pbs_local_num_read_per_job[split_rank] -=  start_index;
		all_reads_offset_source_index = (size_t *)malloc(sizeof(size_t) * total_num_read);

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
		// we gather the offset destination sorted index
		time_count = MPI_Wtime();
		chosen_split_rank_gather_size_t(split_comm,
				split_rank,
				dimensions,
				chosen_split_rank,
				pbs_local_num_read_per_job[split_rank],
				pbs_local_num_read_per_job,
				pbs_start_num_coordinates_per_jobs,
				all_reads_offset_source_index,
				pbs_global_reads_offset_source_index,
				start_index);


		if (split_rank == chosen_split_rank)
			fprintf(stderr,	"rank %d :::::[WRITE] Time to gather all_reads_coordinates_index %f seconds\n", split_rank, MPI_Wtime() - time_count);

		free(pbs_global_reads_offset_source_index); //ok
		free(pbs_local_num_read_per_job);

	} //end if (split_rank < dimensions)

	/*
	 * Bitonic is finished
	 *
	 * Now we re-order and dispatch
	 * the ranks, source offset and sizes
	 *
	 */

	free(pbs_start_num_coordinates_per_jobs); //ok

	// we create a datatype by split_rank

	size_t *all_offset_source_to_send = NULL;

	if (split_rank == chosen_split_rank){

		size_t k;
		//we reorder all_reads_sizes
		all_offset_source_to_send = (size_t *)malloc(sizeof(size_t) * total_num_read);
		all_read_size_to_send = (int *)malloc(sizeof(int) * total_num_read);

		for (k = 0; k < total_num_read ; k++){

			all_offset_source_to_send[k] = all_offset_source_file[all_reads_offset_source_index[k]];
			all_read_size_to_send[k] = all_read_size[all_reads_offset_source_index[k]];

		}

		free(all_offset_source_file);
		free(all_read_size);

	} //end if (split_rank == chosen_rank)

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
	MPI_Barrier(split_comm);
	/*
	 * first we send the offsets of the reads in
	 * the destination file
	 */

	/*
	 * second we send the offsets of the reads in
	 * the source file
	 */

	if (split_rank != master_job){
		//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we recieve offset source \n", rank);
		MPI_Recv(offset_source, local_readNum, MPI_LONG_LONG_INT, master_job, 0, split_comm, &status);

	}
	else {

		size_t ind = start_num_reads_per_jobs[master_job];
		for (k = 0; k < num_reads_per_jobs[master_job]; k++){
			offset_source[k] = all_offset_source_to_send[ind];
			ind++;
		}

		for(j = 0; j < split_size; j++){

			if (j != master_job){
				//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we send all_offset_sourceto_send \n", rank);
				MPI_Send(&all_offset_source_to_send[start_num_reads_per_jobs[j]], num_reads_per_jobs[j], MPI_LONG_LONG_INT, j, 0, split_comm);
			}
		}
	}

	/*
	 * three we send the sizes of the reads in
	 *
	 */

	if (split_rank != master_job){

		//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we recieve new_read_size \n", rank);
		MPI_Recv(size_source, local_readNum, MPI_INT, master_job, 0, split_comm, &status);
	}
	else {

		size_t ind = start_num_reads_per_jobs[master_job];
		for (k = 0; k < num_reads_per_jobs[master_job]; k++){
			size_source[k] = all_read_size_to_send[ind];
			ind++;
		}
		for(j = 0; j < split_size; j++){
			if (j != master_job){
				//fprintf(stderr, "Rank %d ::::: Phase 2 :::: we send all_read_size_to_send \n", rank);
				MPI_Send(&all_read_size_to_send[start_num_reads_per_jobs[j]], num_reads_per_jobs[j], MPI_INT, j, 0, split_comm);
			}
		}
	}


	MPI_Barrier(split_comm);

	if (split_rank == master_job){
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


	//previous version of local sort with output of permutation
	base_arr2 = offset_source;
	//new version of the local sort
	qksort(offset_source_index, local_readNum, sizeof(size_t), 0, local_readNum - 1, compare_size_t);

	/*
	 * reorder the offset and the size
	 * of the reads according to sort
	 */

	for(j = 0; j < local_readNum; j++){
		offset_source_sorted[j] = offset_source[offset_source_index[j]];
		read_size_sorted[j] = size_source[offset_source_index[j]];
	}

	size_t new_data_sz = 0;


	MPI_Barrier(split_comm);

	/*
	 * first we pack data in data_pack1
	 * to have consecutive reads of the
	 * chromosoms
	 */

	int m;
	size_t i;


	MPI_Datatype dt_data0;
	char *data_pack;

	//we compute the size of data_pack
	for (k = 0; k < local_readNum; k++){
		new_data_sz += read_size_sorted[k];
	}

	data_pack = malloc(new_data_sz +1);
	data_pack[new_data_sz] = 0;

	//we create a new offset vector
	//wich is the offset in the data
	size_t *offset_in_data = (size_t *) malloc(local_readNum* sizeof(size_t));

	//we compute the new offset of reads in data buffer
	//we remove the start offset in the file

	for (k=0; k < local_readNum;k++){
		offset_in_data[k] = offset_source[k] - start_offset_in_file;
	}

	MPI_Type_hindexed(local_readNum, &size_source[0], (MPI_Aint*)offset_in_data, MPI_CHAR, &dt_data0);
	MPI_Type_commit(&dt_data0);

	int position=0;

	MPI_Pack(data, 1, dt_data0, data_pack, new_data_sz, &position, COMM_WORLD);

	//Here we are going to send the data to a buffer in the next rank job
	//The next job will compress data and send it back to the prvious job

	//number of block to send
	int blocksize = 1;

	size_t *y_message_sz = (size_t *)calloc(split_size,  sizeof(size_t));
	size_t *y_read_num = (size_t *)calloc(split_size,  sizeof(size_t));

	// copy x into the correct location in y
	y_message_sz[split_rank * blocksize] = new_data_sz;
	y_read_num[split_rank * blocksize] = local_readNum;


	int successor = ( split_rank + 1 ) % split_size;
	int predecessor = ( split_rank - 1 + split_size ) % split_size;


	size_t send_offset;
	size_t recv_offset;

	for (i = 0; i < split_size - 1 ; i++) {

		send_offset = ( ( split_rank - i + split_size ) % split_size );
		recv_offset = ( ( split_rank - i - 1 + split_size ) % split_size );

		// first we size the size of the buffer to send
		MPI_Send( y_message_sz + send_offset, blocksize , MPI_LONG_LONG_INT, successor, 0, split_comm);
		MPI_Recv( y_message_sz + recv_offset, blocksize , MPI_LONG_LONG_INT, predecessor, 0, split_comm, &status);

	}

	//we create a buffer for recieved data
	char *char_buff_uncompressed = malloc((y_message_sz[predecessor] + 1) * sizeof(char));
	char_buff_uncompressed[y_message_sz[predecessor]] = 0;

	//now we send data
	MPI_Sendrecv(data_pack, new_data_sz, MPI_CHAR, successor, 0,
				 char_buff_uncompressed, y_message_sz[predecessor],
				 MPI_CHAR, predecessor, 0, split_comm,  &status);

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
	int compressed_size = 0;

	if (fp->uncompressed_block == NULL)
	   fp->uncompressed_block = malloc(fp->uncompressed_block_size);

	input = (void *)char_buff_uncompressed;
	block_length = fp->uncompressed_block_size;
	bytes_written = 0;
	uint8_t *compressed_buff =  malloc(strlen(char_buff_uncompressed) * sizeof(uint8_t));

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

	free(y_message_sz);

	//we compress the neader


	BGZF *fp_header;
	fp_header = calloc(1, sizeof(BGZF));
	uint8_t *compressed_header = NULL;
	int compressed_size_header = 0;

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

	MPI_Barrier(split_comm);

	//we trade the blocks
	//in this phase the predeccor become the succesor
	size_t compressed_sz_to_send = compressed_size;
	size_t compressed_sz_to_recv = 0;

	// copy x into the correct location in y
	int predecessor_back = ( split_rank + 1 ) % split_size;
	int successor_back = ( split_rank - 1 + split_size ) % split_size;

	// first we size the size of the buffer to send
	MPI_Send( &compressed_sz_to_send, blocksize , MPI_LONG_LONG_INT, successor_back, 0, split_comm);
	MPI_Recv( &compressed_sz_to_recv, blocksize , MPI_LONG_LONG_INT, predecessor_back, 0, split_comm, &status);

	//we create a buffer for recieved data
	uint8_t *buff_compressed = malloc(compressed_sz_to_recv * sizeof(uint8_t));
	//char *buff_compressed = malloc(compressed_sz_to_recv * sizeof(char));
	//now we send data
	MPI_Sendrecv(compressed_buff, compressed_sz_to_send, MPI_UNSIGNED_CHAR, successor_back, 0,
			buff_compressed, compressed_sz_to_recv, MPI_UNSIGNED_CHAR, predecessor_back, 0, split_comm,  &status);

	//fprintf(stderr, "rank %d :::: we recieve from %d a compressed buffer of size %zu \n", rank,
	//		successor, compressed_sz_to_recv );

	MPI_Barrier(split_comm);
	size_t compSize = compressed_sz_to_recv;
	/*
	 * Now we write results of compression
	 */

	MPI_Barrier(split_comm);
	size_t write_offset = 0;

	MPI_Offset * y = (MPI_Offset *) calloc(split_size, sizeof(MPI_Offset));
	MPI_Offset * y2 = (MPI_Offset *) calloc(split_size+1, sizeof(MPI_Offset));

	MPI_Gather(&compSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, split_comm);

	//now we make a cumulative sum
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

	//we do a gather in replacement of the the ring pass
	//fprintf(stderr, "Proc %d:::: we call MPI_Scatter \n", rank);

	MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, split_comm);

	//fprintf(stderr, "Proc %d:::: we recieve offset %zu\n", rank, (size_t)write_offset);

	/*
	 * now we compute offset where to write
	 */

	// we create the path where to write for collective write
	path = (char*)malloc((strlen(output_dir) + strlen(chrName) + 40) * sizeof(char));
	sprintf(path, "%s/%s.bam", output_dir, chrName);


	ierr = MPI_File_open(split_comm, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);

	if (ierr) {
		fprintf(stderr, "Rank %d failed to open %s.\nAborting.\n\n", split_rank, path);
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
	MPI_File_write(out, buff_compressed, (size_t)compSize, MPI_BYTE, &status);
	if (split_rank == master_job)
		fprintf(stderr, "Rank %d :::::[WRITE] Time for chromosome %s writing %f seconds\n", split_rank, chrName, MPI_Wtime()-time_count);

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

	/*
	fprintf(stderr, "Proc %d::::local size to write %zu \n", rank,	dataSize);

	//number of block to send
	MPI_Offset *y = (MPI_Offset *) calloc(num_proc, sizeof(MPI_Offset));
	MPI_Offset *y2 = (MPI_Offset *) calloc(num_proc+1, sizeof(MPI_Offset));
	//we wait all processors
	MPI_Barrier(MPI_split_comm);
	MPI_Gather(&dataSize, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, MPI_split_comm);
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
	MPI_Barrier(MPI_split_comm);
	//we do a gather in replacement of the the ring pass
	MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &write_offset, 1, MPI_LONG_LONG_INT, 0, MPI_split_comm);

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
	ierr = MPI_File_open(MPI_split_comm, path, MPI_MODE_WRONLY  + MPI_MODE_CREATE, finfo, &out);

	if (ierr) {
		fprintf(stderr, "Rank %d failed to open %s.\nAborting.\n\n", rank, path);
		MPI_Abort(split_comm, ierr);
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

	MPI_Barrier(split_comm);
	MPI_File_close(&out);
	//don't close in it will be closed by the discordant part
	//MPI_File_close(&in);
	free(char_buff_uncompressed);
	free(read_size);
	free(read_size_sorted);
	free(offset_source_index);
	free(offset_source_sorted);
	free(offset_source);
	free(size_source);
	free(y);
	free(y2);
	free(path);
	free(data_pack);

}







