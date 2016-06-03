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
     sort_main.c

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/


#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <mpi.h>
#include <assert.h>
#include <fcntl.h>
#include <inttypes.h>
#include <libgen.h>
#include <unistd.h>

#include "mpi_globals.h"
#include "merge.h"
#include "diffuse.h"
#include "preWrite.h"
#include "write2.h"
#include "mergeSort.h"
#include "parser.h"
#include "merge_utils.h"
#include "bufferized_read.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define	PRIoff PRId64

#define	MPI_OFF_T MPI_LONG_LONG_INT

/*
 * Constant for Lustre striping
 * To adapt depending on the file server
 * STRIPING FACTOR is the numer of servers
 * where your files are striped
 *
 * STRIPING UNIT is the size of the stripes
 *
 * if you change those values
 *
 */
#define STRIPING_FACTOR "4"

//#define STRIPING_UNIT "268435456"   // 256 MB
#define STRIPING_UNIT "536870912"   // 500 MB
//#define STRIPING_UNIT "1073741824"  // 1GB
//#define STRIPING_UNIT "1610612736"  // 1.5GB
//#define STRIPING_UNIT "2147483648"  // 2GB
//#define STRIPING_UNIT "2684354560"  // 2.5GB
//#define STRIPING_UNIT "3221225472"  // 3GB
//#define STRIPING_UNIT "3758096384"  // 3.5GB

/*
 * Constant for MPI IO
 * For description of the constant see here
 *
 * https://fs.hlrs.de/projects/craydoc/docs/books/S-2490-40/html-S-2490-40/chapter-sc4rx058-brbethke-paralleliowithmpi.html
 */

#define NB_PROC  "2" //numer of threads for writing
#define CB_NODES "4" //numer of server for writing
#define CB_BLOCK_SIZE  "268435456" /* 256 MBytes - should match FS block size */
#define CB_BUFFER_SIZE  "536870912" /* multiple of the block size by the number of proc*/
#define DATA_SIEVING_READ "enable"

/*
 * Capacity constant no need to change it
 */
#define DEFAULT_MAX_SIZE 6000000000 //Default capacity for one process: 6G
#define DEFAULT_INBUF_SIZE  (1024*1024*1024)


static void usage(const char *);


int main (int argc, char *argv[]){

	DIR *dir = NULL;
	MPI_File mpi_filed;

	MPI_Offset unmapped_start, discordant_start;
	int num_proc, rank;
	int nbchr, i, paired=0; //we assume the reads are single ended
	int ierr, errorcode = MPI_ERR_OTHER;
	char *file_name, *output_dir, *sam_dir;

	char *header;
	char **chrNames;
	struct stat fileSize;
	unsigned int headerSize;
	unsigned char threshold = 0;
	char sender;

	size_t unmappedSize = 0;
	size_t array_max_size = DEFAULT_MAX_SIZE; //maximum amount of data sent to the father
	size_t *readNumberByChr = NULL, *localReadNumberByChr = NULL;
	Read **reads;
	size_t *count_diffuse = NULL;
	size_t **send_diffuse;
	size_t *dsend[3];

	//init dsend


	clock_t tic, toc;
	int compression_level = 3; //by default compression is 3

	char *rbuf;
	size_t fsiz, lsiz, loff, *goff;

	MPI_Info finfo;
	MPI_Init(&argc,&argv);
	/*
	if (argc < 4){
		fprintf(stderr, "Invalid arguments.\nShutting down.\n");
		MPI_Finalize();
		return 0;
	}
	*/
    //finds out process rank
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//finds out number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	if(!rank)fprintf(stderr,"Number of processes : %d\n",num_proc);

	file_name = argv[1];
	output_dir = argv[2];

	if(rank == 0)fprintf(stderr, "file to read : %s\n", file_name);
	if(rank == 0)fprintf(stderr, "output : %s\n", output_dir);

	if (argc < 2) {
		usage(basename(*argv));
		MPI_Finalize();
		return 1;

	}

	//looking for optionsinfo
	for(i = 0; i < argc; i++){
		if(argv[i][0] == '-'){
			if(argv[i][1] == 'q'){
				threshold = atoi(argv[i+1]);
				if(!rank)fprintf(stderr, "Reads' quality threshold : %d\n", threshold);
			}

			if(argv[i][1] == 'p'){
				paired = 1;
				if(!rank)fprintf(stderr, "Reads are paired \n");
			}

			if(argv[i][1] == 'c'){
				compression_level = atoi(argv[i+1]);
				if(!rank)fprintf(stderr, "Compression Level is : %d\n", compression_level);
			}

			if(argv[i][1] == 'h'){
				usage(basename(*argv));
				MPI_Finalize();
				return 0;
			}

		}
	}

	//checking if everything's fine with the output directory, shutting down otherwise
	sam_dir = (char*)malloc(strlen(output_dir)+11);
	sprintf(sam_dir, "%s", output_dir);
	if(!rank){
		dir = opendir(sam_dir);
		if(!dir){
			perror("Failed to open directory");
			if(errno == ENOENT){
				fprintf(stderr, "rank %d Making directory...\n", rank);
				ierr = mkdir(sam_dir, S_IRWXU);
				if(ierr){
					perror("Failed to make directory");
					fprintf(stderr, "Shutting down...\n");
					MPI_Abort(MPI_COMM_WORLD, errorcode);
					exit(2);
				}
			}
			else{
				fprintf(stderr, "Shutting down...\n");
				MPI_Abort(MPI_COMM_WORLD, errorcode);
				exit(2);
			}
		}
		else{
			closedir(dir);
		}
	}


	//task FIRST FINE TUNING FINFO FOR READING OPERATIONS

	/*
	 * In this part you shall adjust the striping factor and unit according
	 * to the underlying filesystem
	 */
	MPI_Info_create(&finfo);
	MPI_Info_set(finfo,"striping_factor", STRIPING_FACTOR);
	MPI_Info_set(finfo,"striping_unit", STRIPING_UNIT); //2G striping
	MPI_Info_set(finfo,"ind_rd_buffer_size", STRIPING_UNIT); //2gb buffer
	MPI_Info_set(finfo,"romio_ds_read",DATA_SIEVING_READ);
		
	/*
	 * for collective reading and writing
	 * should be adapted too and tested according to the file system
	 */
	MPI_Info_set(finfo,"nb_proc", NB_PROC);
	MPI_Info_set(finfo,"cb_nodes", CB_NODES);
	MPI_Info_set(finfo,"cb_block_size", CB_BLOCK_SIZE);
	MPI_Info_set(finfo,"cb_buffer_size", CB_BUFFER_SIZE);
	

	//we open the input file
	ierr = MPI_File_open(MPI_COMM_WORLD, file_name,  MPI_MODE_RDONLY , finfo, &mpi_filed);
	//assert(in != -1);
	if (ierr){
		if (rank == 0) DEBUG("%s: Failed to open file in process 0 %s\n", argv[0], argv[1]);
		MPI_Abort(MPI_COMM_WORLD, errorcode);
		exit(2);
	}

	//then we get the file size
	size_t input_file_size = stat(file_name, &fileSize);

	if (input_file_size == -1){
		fprintf(stderr,"Failed to find file size\n");
		MPI_Abort(MPI_COMM_WORLD, errorcode);
		exit(0);
	}

	input_file_size = (long long)fileSize.st_size;
	if(!rank)fprintf(stderr, "The size of the file is %zu\n", fileSize.st_size);

	/* Get chunk offset and size */
	fsiz = fileSize.st_size;
	lsiz = fsiz / num_proc;
	loff = rank * lsiz;
	size_t lsiz2 = 150*sizeof(char)*1000; // load only the begining of file to check the read name and header
	if(lsiz2>lsiz){
		lsiz2 = lsiz;
	}

	rbuf = (char*)malloc((lsiz2 + 1)*sizeof(char));
	tic = MPI_Wtime();
	MPI_File_read_at(mpi_filed, loff, rbuf, lsiz2, MPI_CHAR, MPI_STATUS_IGNORE);
	fprintf(stderr, "%d (%.2lf)::::: ***GET OFFSET BLOCKS TO READ ***\n", rank,(double)(MPI_Wtime()-tic));


	//FIND HEADERSIZE AND CHRNAMES AND NBCHR
	tic = MPI_Wtime();
	headerSize=find_header(rbuf,rank,&unmappedSize,&nbchr,&header,&chrNames);
	fprintf(stderr, "%d (%.2lf)::::: ***HEADER PART%d***\n", rank,(double)(MPI_Wtime()-tic),headerSize);
	free(rbuf);

	//We place file offset of each process to the begining of one read's line
	goff=init_goff(mpi_filed,headerSize,fileSize.st_size,num_proc,rank);

	//We calculate the size to read for each process
	lsiz = goff[rank+1]-goff[rank];
	//NOW WE WILL PARSE
	int j=0;
	size_t poffset = goff[rank]; //Current offset in file sam

	//nbchr because we add the discordant reads in the structure
	reads = (Read**)malloc((nbchr)*sizeof(Read));//We allocate a linked list of struct for each Chromosome (last chr = unmapped reads)
	readNumberByChr = (size_t*)malloc((nbchr)*sizeof(size_t));//Array with the number of reads found in each chromosome
	localReadNumberByChr = (size_t*)malloc((nbchr)*sizeof(size_t));//Array with the number of reads found in each chromosome
	Read ** anchor = (Read**)malloc((nbchr)*sizeof(Read));//Pointer on the first read of each chromosome

	//Init first read
	for(i = 0; i < (nbchr); i++){
		reads[i] = malloc(sizeof(Read));
		reads[i]->coord = 0;
		anchor[i] = reads[i];
		readNumberByChr[i]=0;
	}

	toc = MPI_Wtime();
	char *local_data = NULL; //Where we load file sam

	//We read the file sam and parse
	while(poffset < goff[rank+1]){

		size_t size_to_read = 0;

		//Reading in multiple times because of MPI_File_read limits
		if( (goff[rank+1]-poffset) < DEFAULT_INBUF_SIZE ){
			size_to_read = goff[rank+1]-poffset;
		}
		else{
			size_to_read = DEFAULT_INBUF_SIZE;
		}

		// we load the buffer
		local_data=(char*)calloc(size_to_read+1,sizeof(char));

		/*
		 * TODO: Issue with MPI_BOTTOM on certain infrastructure problem
		 * 		 of randomization with (void *)0
		 * 		 We comme back to the previous mpi_file_read_at line 311
		 *
		 */

		// Original reading part is before 18/09/2015
		MPI_File_read_at(mpi_filed, (MPI_Offset)poffset, local_data, size_to_read, MPI_CHAR, MPI_STATUS_IGNORE);
		size_t local_offset=0;


		// modification 18/09/2015
		// we create a Datatype for the block
		// creation of datatype
		// Variable for datatype struct almost classic
		// The idea is to decompose the block to read
		// in small chunks (of integer size)
		// for 1gb there's 1.5 millions read
		// creation of datatype

		/*
		MPI_Datatype dt_view;
		MPI_Datatype dt_data;

		size_t num_data_type_block =150;
	 	int k=0;
	 	int blocklens[num_data_type_block];
	 	MPI_Aint indices[num_data_type_block];
	 	MPI_Datatype oldtypes[num_data_type_block];

	 	size_t *vector_offset = (size_t*)malloc(num_data_type_block*sizeof(size_t));

	 	//We divide the size_to_read in readNum small size
	 	int chunck_size = size_to_read/num_data_type_block;
	 	int rest = size_to_read - (num_data_type_block* chunck_size);

	 	//Allocate data
	 	for(k = 0; k < size_to_read; k++){
	 		local_data[k] = 0;
	 	}
	 	char *u = local_data;
	 	for ( k = 0; k < num_data_type_block; k++){
	 		blocklens[k] = chunck_size;

	 		//MPI_Get_address((u + k*chunck_size), &indices[k]);
	 		indices[k] = (MPI_Aint)(u + k*chunck_size);
	 		vector_offset[k] = poffset + k*chunck_size;
	 		oldtypes[k] = MPI_CHAR;
	 	}
	 	blocklens[num_data_type_block -1]+=rest;

	 	size_t res=0;
	 	for ( k = 0; k < num_data_type_block; k++){
	 		res+=blocklens[k];
	 	}
	 	assert(res==size_to_read);

	 	MPI_Type_create_struct(num_data_type_block, blocklens, indices, oldtypes, &dt_data);
	 	MPI_Type_commit(&dt_data);

		MPI_Type_create_hindexed(1, blocklens, (MPI_Aint *)vector_offset, MPI_CHAR, &dt_view);
		MPI_Type_commit(&dt_view);
		//TODO: see if initialization is needed
		MPI_File_set_view(mpi_filed, 0, MPI_CHAR, dt_view, "native", finfo);
		//MPI_File_read(mpi_filed, local_data, 1, MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_read(mpi_filed, MPI_BOTTOM, 1, dt_data, MPI_STATUS_IGNORE);

		MPI_Type_free(&dt_data);
		MPI_Type_free(&dt_view);
		 */

		////////////////////////////////////
		///// fin modification 18/09/2015
		////////////////////////////////////


		//we look where is the last line read for updating next poffset
		size_t offset_last_line = size_to_read-1;
		while(local_data[offset_last_line] != '\n'){
			offset_last_line -- ;
		}
		//If it s the last line of file, we place a last '\n' for the function tokenizer
		if(rank == num_proc-1 && ((poffset+size_to_read) == goff[num_proc])){
			local_data[offset_last_line]='\n';
		}

		//Now we parse Read in local_data
		parser_paired(local_data, rank, poffset, threshold, nbchr, &readNumberByChr, chrNames, &reads);

		//we go to the next line
		poffset+=(offset_last_line+1);
		local_offset+=(offset_last_line+1);
		free(local_data);
	}

	//MPI_File_seek(mpi_filed, goff[rank], MPI_SEEK_SET);

	fprintf(stderr, "%d (%.2lf)::::: *** FINISH PARSING FILE ***\n", rank,(double)(MPI_Wtime()-toc));

	free(goff);
	//MPI_File_close(&mpi_filed);
	MPI_Barrier(MPI_COMM_WORLD);

	//We set attribute next of the last read and go back to first read of each chromosome
	for(i = 0; i < nbchr; i++){
		reads[i]->next = NULL;
		reads[i] = anchor[i];
	}
	free(anchor);

	//We count how many reads we found
	size_t nb_reads_total =0,nb_reads_global =0;
	for(j=0;j<nbchr;j++){
		nb_reads_total+=readNumberByChr[j];
	}

	MPI_Allreduce(&nb_reads_total, &nb_reads_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	fprintf(stderr, "Number of reads on rank %d = %zu/%zu \n", rank, nb_reads_total, nb_reads_global);


	/*
	 *
	 * In this part we are going to define a new communicator
	 * in the case we have a job without reads to sort.
	 *
	 * We chose to use communicator to avoid modification
	 * in the communication part of the writing process that could be
	 * exhausting
	 *
	 * To do that:
	 *
	 * 1) we rank 0 collects all the local number of reads to sort.
	 * 2) the rank 0 compute the new rank in the new communicator
	 * 3) the rank 0 broadcast the new rank rank
	 * 4) MPI_COMM_SPLIT create the new communicator NEW_MPI_COMM_WORLD
	 * 5) free the new communicator
	 * 6) loop the next chromosome
	 *
	 *
	 */
	/*
	int chr_num = 0;

	for (chr_num = 0; chr_num < (nbchr-2); chr_num++){

		int *localReadsNum_vec = (int *) calloc(num_proc);
		int *jobRank_vec = (int *) calloc(num_proc); //vector of rank jobs before split comm

		//we build a vector with rank job
		int zero_counter_1 = 0; // counter on localReadsNum_vec
		int zero_counter_2 = 0; // counter on newRank_vec
		int j = 0;

		//rank 0 gather the vector
		MPI_Gather( localReadNumberByChr[chr_num], 1, MPI_LONG_LONG_INT, localReadsNum_vec[chr_num], 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD)

		if (rank == 0) {
			//we compute the number of 0 in the localReadsNum_vec
			for(i = 0; i < num_proc; i++){
				if (localReadsNum_vec[i] == 0) counter++;
				jobRank_vec[i] = i;
			}
			// if no jobs without reads we do nothing
			if ( zero_counter_1 == 0 ){
				// nothing to do
				NEW_MPI_COMM_WORLD = MPI_COMM_WORLD;
				break; //we exit the for loop
			}
			else{

				int jobRank_vec_sz = num_proc - zero_counter_1;
				// allocate new rank vector size of counter
				int *jobRank_vec = (int *) calloc(jobRank_vec_sz);
				// we initialize jobRank_vec
				for(i = 0; i < jobRank_vec_sz; i++){
					jobRank_vec[i] = i;
				}

				// if the localnumRead is 0 we update the rank
				for (j = 0; j < zero_counter_1; j++){

					if ( localReadsNum_vec[j] == 0 ){

					}
					else{
					}

				}

		int key = 0; //everybody has a zero key
		int color = 0;

		//first we mus

		if (localReadNumberByChr[i] == 0){
			key=rank; 				 //except the rank who has no reads to sort
			color=MPI_UNDEFINED;     //we exclude this rank
		}

		int color = 0;
		MPI_Comm* NEW_MPI_COMM_WORLD;

		if (localReadNumberByChr[i] == 0){
			// then we shall split the communicator to exclude the rank with
			// no reads
			MPI_Comm_split( MPI_COMM_WORLD, color, key, NEW_MPI_COMM_WORLD);

			//now we change the rank in the reads structure

			//we must update the num_proc value
			MPI_Comm_rank(New_Comm, &new_id);
			//we must update the rank

		}

	}
	*/


	/*
	*  We write the mapped reads in a file named chrX.bam
	*
	*/

	for(i = 0; i < (nbchr-2); i++){

		localReadNumberByChr[i] = readNumberByChr[i];

		if(reads[i] && reads[i]->next && reads[i]->next->next){
			mergeSort(reads[i], readNumberByChr[i]);
		}

		reads[i] = reads[i]->next;
		indexing(rank, readNumberByChr[i], reads[i], dsend);

		count_diffuse = NULL; //used for the formatting
		sender = merge(rank, num_proc, headerSize, readNumberByChr[i], array_max_size, &count_diffuse, &send_diffuse, dsend);

		size_t *offsets = (size_t*)malloc(localReadNumberByChr[i]*sizeof(size_t));

		diffuse(offsets, rank, num_proc, sender, localReadNumberByChr[i], count_diffuse, send_diffuse);

		writeSam(rank, output_dir, header, localReadNumberByChr[i], chrNames[i],
				reads[i], offsets, num_proc, MPI_COMM_WORLD, file_name, mpi_filed, finfo, compression_level);

		MPI_Barrier(MPI_COMM_WORLD);
	}


	/*
	 *  We write the unmapped reads in a file named unmapped.bam
	 *
	 */

	reads[nbchr-2] = reads[nbchr-2]->next;
	localReadNumberByChr[nbchr-2] = readNumberByChr[nbchr-2];
	unmapped_start = unmappedOffset(rank, num_proc, unmappedSize, headerSize, nbchr-2, localReadNumberByChr[nbchr-2]);

	if(!unmapped_start){
		fprintf(stderr, "No header was defined.\nShutting down.\n");
		MPI_Finalize();
		return 0;
	}
	else{

		writeSam_unmapped(rank, output_dir, header, localReadNumberByChr[nbchr-2], chrNames[nbchr-2],
			reads[nbchr-2], num_proc, MPI_COMM_WORLD, file_name, mpi_filed, finfo, compression_level);

		while( reads[nbchr-2]->next != NULL){
			Read *tmp_chr = reads[nbchr-2];
			reads[nbchr-2] = reads[nbchr-2]->next;
			free(tmp_chr->next);
		}
	}

	/*
     *  Now we write the discordant reads in a file names discordant_reads.sam
	 */

	reads[nbchr-1] = reads[nbchr-1]->next;
	localReadNumberByChr[nbchr-1] = readNumberByChr[nbchr-1];
	discordant_start = discordantOffset(rank, num_proc, unmappedSize, headerSize, nbchr-1, localReadNumberByChr[nbchr-1]);

	if(!discordant_start){
		fprintf(stderr, "No header was defined.\nShutting down.\n");
		MPI_Finalize();
		return 0;
	}
	else{

		/*
		 * discordant reads are reads where pair read don't align
		 * in the same chromosome or one read is aligned and not the pair
		 */

		writeSam_discordant(rank, output_dir, header, localReadNumberByChr[nbchr-1], chrNames[nbchr-1],
			reads[nbchr-1], num_proc, MPI_COMM_WORLD, file_name, mpi_filed, finfo, compression_level);

		while( reads[nbchr-1]->next != NULL){
			Read *tmp_chr = reads[nbchr-1];
			reads[nbchr-1] = reads[nbchr-1]->next;
			free(tmp_chr->next);
		}
		//free(reads[i]);
		if(!rank){
			free(header);
		}
		fprintf(stderr,"Rank %d finished.\n", rank);
		printf("rank %d \n", rank);
		MPI_Finalize();

	}


	free(localReadNumberByChr);

	for(i = 0; i < nbchr; i++)
		free(chrNames[i]);
	free(chrNames);
	// task: FREE READS
	//free(reads);
	free(readNumberByChr);
	return 0;
}


void create_read_dt_for_parser(int rank, int num_proc, int *ranks, int* buffs, char** data, MPI_Datatype* dt, size_t readNum)
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
 	MPI_Datatype oldtypes[readNum];

 	/* Adress originally referencing on data
 	 * data : char** buff in which the read data must be organized by destination rank
 	 */
 	MPI_Aint adress_to_write_in_data_by_element[num_proc];

 	//init adress_to_write_in_data_by_element
 	for(i = 0; i < num_proc; i++){
 		adress_to_write_in_data_by_element[(rank+i)%num_proc] = (MPI_Aint)data[(rank-i+num_proc)%num_proc];
 	}

 	//double time_count = MPI_Wtime();
 	//Set all classic datatype values
 	for(i = 0; i < readNum; i++){
 		/*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TRICKY PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 		 * Basically what we say here is that the i(th) item that we are going to read goes to the data row corresponding to it's destination rank
 		 *
 		 * indices[i] : The adress for the datatype to which the i(th) item should be written.
 		 * 			Generally this is a relative adress. Like 0, 8, 16.
 		 * 			Here, we use the adress in memory in order to be able to use an array of array without any problem.
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
 	//time_count = MPI_Wtime();
 	MPI_Type_create_struct(readNum, blocklens, indices, oldtypes, dt);
 	//fprintf(stderr, "Rank %d :::::[create_read_dt] Time for creating struct %f seconds\n", rank, MPI_Wtime()-time_count);
 }


static void usage(const char *prg) {

	fprintf(stderr, "Program: MPI version for sorting FASTQ data\n"
		"Version: v1.0\n"
		"Contact 1: Frederic Jarlier (frederic.jarlier@curie.fr) \n"
		"usage : mpirun -n TOTAL_PROC %s FILE_TO_SORT OUTPUT_FILE -q QUALITY \n"
		"output : a bam files per chromosome, a bam file of unmapped reads \n"
		"                 a bam files of discordants reads. \n"
		"Discordants reads are reads where one pairs align on a chromosome \n"
		"and the other pair align on another chromosome \n"
		"Unmapped reads are reads without coordinates on any genome \n"
		"Requirements :  For perfomances matters the file you want to sort has to be \n"
		"       stripped on parallel file system. \n"
		"With Lustre you mention it with lfs setstripe command like this \n"
		"lfs set stripe -c stripe_number -s stripe_size folder           \n"
		"In mpiSort you mention the number of stripes with -c options \n"
		"and mention the sripe size with with -s options. Example of command     \n"
		"If you file is striped with 16 servers and chunk size of 1Gb the lfs options will be -c 16 -s 1 \n",
		prg);

	return; }
