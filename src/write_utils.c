#include "write_utils.h"

size_t g_wu_readNum;
size_t g_wu_totalReadNum;
size_t g_wu_master;
MPI_Comm COMM_WORLD;
MPI_Status status;

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

 	//fprintf(stderr, "Rank %d :::::[WRITE] readNum = %zu \n", rank, readNum );
 	MPI_Barrier(COMM_WORLD);
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


void send_size_t_all_to_master(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs,
		size_t *start_size_per_job, size_t *all_data, size_t *data)
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
				MPI_Recv(temp_buf, size_per_jobs[j], MPI_LONG_LONG_INT, j, 0, COMM_WORLD, &status);

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

		MPI_Send(data, size, MPI_LONG_LONG_INT, master,  0, COMM_WORLD);
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
				MPI_Recv(temp_buf, size_per_jobs[j], MPI_LONG_LONG_INT, j, 0, COMM_WORLD, &status);

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

		MPI_Send(data, size, MPI_LONG_LONG_INT, master,  0, COMM_WORLD);
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
				MPI_Recv(temp_buf, size_per_jobs[j], MPI_LONG_LONG_INT, j, 0, COMM_WORLD, &status);
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

		MPI_Send(data + start_index, size, MPI_LONG_LONG_INT, master,  0, COMM_WORLD);
	}
}



void send_size_t_master_to_all(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs,
		size_t *start_size_per_job, size_t *all_data, size_t *data)
{
	MPI_Status status;
	int j, k;
	if (rank != master){
		//fprintf(stderr, "%d ::::: [send_size_t_master_to_all] rank %d recv %zu from %d \n",rank, rank, size, master);

		MPI_Recv(data, size, MPI_LONG_LONG_INT, master, 0, COMM_WORLD, &status);


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
						size_per_jobs[j], MPI_LONG_LONG_INT, j, 0, COMM_WORLD);

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
				MPI_Recv(temp_buf, size_per_jobs[j], MPI_INT, j, 0, COMM_WORLD, &status);

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

		MPI_Send(data, size, MPI_INT, master,  0, COMM_WORLD);
	}
}

void send_int_master_to_all(int rank, int num_proc, int master, size_t size, size_t *size_per_jobs, size_t *start_size_per_job, int *all_data, int *data)
{
	MPI_Status status;
	int j, k;

	if (rank != master){
		MPI_Recv(data, size, MPI_INT, master, 0, COMM_WORLD, &status);
	}
	else {

		size_t ind = start_size_per_job[master];
		for (k = 0; k < size_per_jobs[master]; k++){
			data[k] = all_data[ind];
			ind++;
		}

		for(j = 0; j < num_proc; j++){

			if (j != master){

				//fprintf(stderr, "%d ::::: [send_int_master_to_all] rank %d send %zu to %d from %zu\n",
				//					rank, rank, size_per_jobs[j], j, start_size_per_job[j]);
				MPI_Send(&all_data[start_size_per_job[j]], size_per_jobs[j], MPI_INT, j, 0, COMM_WORLD);
			}
		}
	}
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



void assert_sorted_size_t(size_t * buffer, size_t number)
{
	size_t i;
	for (i = 0; i < number; i++) assert( buffer[i] < buffer[i + 1]);
}
