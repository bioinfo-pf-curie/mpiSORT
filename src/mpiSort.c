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

#define _GNU_SOURCE

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

#include <mpi.h>

#include "diffuse.h"
#include "mergeSort.h"
#include "merge_utils.h"
#include "parser.h"
#include "preWrite.h"
#include "write.h"
#include "mpiSort_utils.h"
#include "write_utils.h"
#include "qksort.h"

#if 0
#include "mpi_globals.h"
#include "merge.h"
#include "diffuse.h"
#include "preWrite.h"
#include "write2.h"
#include "mergeSort.h"
#include "parser.h"
#include "merge_utils.h"
#include "bufferized_read.h"
#endif


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

#define NB_PROC  "40" //numer of threads for writing
#define CB_NODES "4" //numer of server for writing
#define CB_BLOCK_SIZE  "268435456" /* 256 MBytes - should match FS block size */
#define CB_BUFFER_SIZE  "536870912" /* multiple of the block size by the number of proc*/
#define DATA_SIEVING_READ "enable"

/*
 * Capacity constant no need to change it
 */
#define DEFAULT_MAX_SIZE 6000000000 //Default capacity for one process: 6G
#define DEFAULT_INBUF_SIZE  (1024*1024*1024)

/* Maximum chromosome number */
#define MAXNBCHR 256

static void usage(const char *);



void bruck_offsets_dest(int rank, int num_proc, int local_readNum,
		size_t* number_of_reads_by_procs, size_t ** data_offsets,
			int *new_rank, size_t* new_offset, MPI_Comm split_comm)
   {
   	MPI_Comm comm = split_comm;

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

   				//data_offsets[recv_index[m]] = realloc(data_offsets[recv_index[m]], sizeof(size_t)*(recv_size_by_proc[m]));
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



int main (int argc, char *argv[]){

	char *x, *y, *z, *xbuf, *hbuf, *chrNames[MAXNBCHR];
	int fd;
	off_t hsiz;
	struct stat st;

	MPI_File mpi_filed;
	MPI_File mpi_file_split_comm;

	MPI_Offset fileSize, unmapped_start, discordant_start;
	int num_proc, rank;
	int res, nbchr, i, paired=0; //we assume the reads are single ended
	int ierr, errorcode = MPI_ERR_OTHER;
	char *file_name, *output_dir;

	char *header;

	unsigned int headerSize;
	unsigned char threshold = 0;

	size_t input_file_size;
	size_t unmappedSize = 0;
	size_t discordantSize = 0;
	size_t *readNumberByChr = NULL, *localReadNumberByChr = NULL;
	Read **reads;

	double time_count;
	double time_count1;

	MPI_Comm split_comm; //used to split communication when jobs have no reads to sort
	int split_rank, split_size; //after split communication we update the rank

	clock_t tic, toc;
	int compression_level = 3; //by default compression is 3

	char *rbuf;
	size_t fsiz, lsiz, loff, *goff;

	MPI_Info finfo;

	parse_mode = MODE_OFFSET;

	/* Check command line */
	while ((i = getopt(argc, argv, "nc:hpq:")) != -1) {
		switch(i) {
			case 'c': /* Compression level */
				compression_level = atoi(optarg);
				break;
			case 'h': /* Usage display */
				usage(basename(*argv));
				return 0;
			case 'p': /* Paired reads */
				paired = 1;
				break;
			case 'q': /* Quality threshold */
				threshold = atoi(optarg);
				break;
			case 'n':
				parse_mode = MODE_NAME;
				break;
			default:
				usage(basename(*argv));
				return 1;
		}
	}
	if (argc - optind != 2) {
		usage(basename(*argv));
		return 1;
	}
	file_name = argv[optind];
	output_dir = argv[optind+1];

	/* Check arguments */
	res = access(file_name, F_OK|R_OK);
	if (res == -1)
		err(1, "%s", file_name);
	res = access(output_dir, F_OK|W_OK);
	if (res == -1)
		err(1, "%s", output_dir);

	/* MPI inits */
	res = MPI_Init(&argc, &argv);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	assert(res == MPI_SUCCESS);

	g_rank = rank;
	g_size = num_proc;

	/* Small summary */
	if (rank == 0) {
		fprintf(stderr, "Number of processes : %d\n", num_proc);
		fprintf(stderr, "Reads' quality threshold : %d\n", threshold);
		fprintf(stderr, "Compression Level is : %d\n", compression_level);
		fprintf(stderr, "SAM file to read : %s\n", file_name);
		fprintf(stderr, "Output directory : %s\n", output_dir);
	}

	/* Process input file */
	fd = open(file_name, O_RDONLY, 0666);
	assert(fd != -1);
	assert(fstat(fd, &st) != -1);
	xbuf = mmap(NULL, st.st_size, PROT_READ, MAP_FILE|MAP_PRIVATE, fd, 0);
	assert(xbuf != MAP_FAILED);

	/* Parse SAM header */
	memset(chrNames, 0, sizeof(chrNames));
	x = xbuf; nbchr = 0;
	while (*x == '@') {
		y = strchr(x, '\n');
		z = x; x = y + 1;
		if (strncmp(z, "@SQ", 3) != 0) continue;
		/* Save reference names */
		y = strstr(z, "SN:");
		assert(y != NULL);
		z = y + 3;
		while (*z && !isspace((unsigned char)*z)) z++;
		chrNames[nbchr++] = strndup(y + 3, z - y - 3);
		assert(nbchr < MAXNBCHR - 2);
	}
	chrNames[nbchr++] = strdup(UNMAPPED);
	chrNames[nbchr++] = strdup(DISCORDANT);
	hsiz = x - xbuf; hbuf = strndup(xbuf, hsiz);
	fprintf(stderr, "Header has %d+2 references\n", nbchr - 2);

	assert(munmap(xbuf, st.st_size) != -1);
	assert(close(fd) != -1);

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
		if (rank == 0) fprintf(stderr, "%s: Failed to open file in process 0 %s\n", argv[0], argv[1]);
		MPI_Abort(MPI_COMM_WORLD, errorcode);
		exit(2);
	}
	ierr = MPI_File_get_size(mpi_filed, &fileSize);
	assert(ierr == MPI_SUCCESS);
	input_file_size = (long long)fileSize;
	if (rank == 0)
		fprintf(stderr, "The size of the file is %zu\n", input_file_size);

	/* Get chunk offset and size */
	fsiz = input_file_size;
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

	if (parse_mode == MODE_OFFSET)
		asprintf(&header, "@HD\tVN:1.0\tSO:coordinate\n%s", hbuf);
	else {

		header = (char*)malloc(strlen(hbuf)+1);
		memcpy(header, hbuf, strlen(hbuf));
		header[strlen(hbuf)] = 0;

	}

	headerSize = unmappedSize = discordantSize = strlen(header);
	free(hbuf); free(rbuf);

	//We place file offset of each process to the begining of one read's line
	goff=init_goff(mpi_filed,headerSize,input_file_size,num_proc,rank);

	//We calculate the size to read for each process
	lsiz = goff[rank+1]-goff[rank];
	//NOW WE WILL PARSE
	size_t j=0;
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
		// Variable for datatype classic struct
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

	fprintf(stderr, "%d (%.2lf)::::: *** FINISH PARSING FILE chr1:%zu ***\n", rank,(double)(MPI_Wtime()-toc), readNumberByChr[0]);

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
	 *  We write the mapped reads in a file named chrX.bam
	 *
	 */
	MPI_Barrier(MPI_COMM_WORLD);

	for(i = 0; i < (nbchr-2); i++){

		int i1,i2;
		//size_t *localReadsNum_rank0 = (size_t *)malloc(num_proc*sizeof(size_t));
		size_t localReadsNum_rank0[num_proc];
		int file_pointer_to_free = 0;
		int split_comm_to_free = 0;
		//we build a vector with rank job
		int val_tmp1 = 0;
		int val_tmp2 = 0;
		int chosen_rank = 0; //needed to tell what rank is going to compute the color and key
		int chosen_split_rank= 0; //the rank that collect data once the communication splitted normally this rank is 0

		// the color tells in what communicator the rank pertain
		// color = 0 will be the new communicator color
		// otherwise the color is 1
		int *color_vec_to_send =  (int *)malloc(num_proc*sizeof(int));
		// the key value tell the order in the new communicator
		int *key_vec_to_send =  (int *)malloc(num_proc*sizeof(int));

		// first we test if the there's reads to sort
		// rank 0 recieve the sum of all the reads count
		size_t total_reads_by_chr = 0;
		MPI_Allreduce(&readNumberByChr[i], &total_reads_by_chr, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

		//if (total_reads_by_chr == 0)
		//	continue; //pass to next chromosome

		//rank 0 gather the vector
		MPI_Allgather(&readNumberByChr[i] , 1, MPI_LONG_LONG_INT, localReadsNum_rank0 , 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);

		if (rank == 0){
			//we must chose the first rank with reads to sort
			i1=0;
			while ((localReadsNum_rank0[i1] == 0) && (i1 < num_proc)){
				chosen_rank++;
				i1++;
			}
			fprintf(stderr, "rank %d :::: chosen rank = %d \n", rank, chosen_rank);
		}

		//we broadcast the chosen rank
		//task: replace the broadcast with a sendrecieve
		MPI_Bcast( &chosen_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		if (((rank == chosen_rank) || rank == 0) && (chosen_rank != 0)){

			//first we exchange the size o
			if (rank == chosen_rank){
				header=(char *)malloc((headerSize + 1)*sizeof(char));
				MPI_Recv(header, headerSize + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			if (rank == 0){
				MPI_Send(header, headerSize + 1, MPI_CHAR, chosen_rank,  0, MPI_COMM_WORLD);
			}
		}
		else {
			//we do nothing here
			//fprintf(stderr, "rank %d :::: nothing here \n", rank );
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (rank == chosen_rank) {
			int counter = 0;
			//we compute the number of 0 in the localReadsNum_vec
			for(i1 = 0; i1 < num_proc; i1++){
				if (localReadsNum_rank0[i1] == 0) {
						counter++;
					}
			}
			// if no jobs without reads we do nothing
			if ( counter == 0 ){
				// nothing to do we associate split_comm with
				fprintf(stderr, "rank %d :::: we don't split the rank \n", rank);
				split_comm = MPI_COMM_WORLD;
				for (i2 = 0; i2 < num_proc; i2++) {
					if (localReadsNum_rank0[i] == 0) {
						color_vec_to_send[i2] = 1;
						key_vec_to_send[i2] = val_tmp2;
						val_tmp2++;
					} else {
						color_vec_to_send[i2] = 0;
						key_vec_to_send[i2] = val_tmp1;
						val_tmp1++;
					}
				}
			}
			else{
				// now we compute the color according to
				// the number of reads to sort
				fprintf(stderr, "rank %d :::: we split the rank \n", rank);
				for(i2 = 0; i2 < num_proc; i2++){
					if (localReadsNum_rank0[i2] == 0){
						color_vec_to_send[i2] = 1;
						key_vec_to_send[i2] = val_tmp2;
						val_tmp2++;
					} else{
						color_vec_to_send[i2] = 0;
						key_vec_to_send[i2] = val_tmp1;
						val_tmp1++;
					}
				} // end for loop
			}// end if
		}// end if (rank == plit_rank)
		MPI_Barrier(MPI_COMM_WORLD);
		//we create key and color variable for each job
		int local_color = 0;
		int local_key = 0;

		// rank 0 scatter the color and the key vector
		MPI_Scatter( color_vec_to_send, 1, MPI_INT, &local_color, 1, MPI_INT, chosen_rank, MPI_COMM_WORLD);
		MPI_Scatter( key_vec_to_send, 1, MPI_INT, &local_key, 1, MPI_INT, chosen_rank, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		// now we create a communicator
		// we group all communicator
		// with color of zero
		if (local_color == 0){
			MPI_Comm_split( MPI_COMM_WORLD, local_color, local_key, &split_comm);
			ierr = MPI_File_open(split_comm, file_name,  MPI_MODE_RDONLY, finfo, &mpi_file_split_comm);
			//we ask to liberate file pointer
			file_pointer_to_free = 1;
			//we ask to liberate the split_comm
			split_comm_to_free = 1;
		}
		else{
			MPI_Comm_split( MPI_COMM_WORLD, MPI_UNDEFINED, local_key, &split_comm);
			mpi_file_split_comm = mpi_filed;
		}

		//now we change the rank in the reads structure
		if (local_color == 0){

			MPI_Comm_rank(split_comm, &split_rank);
			MPI_Comm_size(split_comm, &split_size);

			//we update g_rank
			g_rank = split_rank;
			g_size = split_size;
		}
		else{
			g_rank = split_rank = rank;
			g_size = split_size = num_proc;
		}


		localReadNumberByChr[i] = readNumberByChr[i];

		MPI_Barrier(MPI_COMM_WORLD);

		if ((local_color == 0) && (i < (nbchr - 2))) {

			//we care for chromosoms reads

			//we do a local merge sort
			if(reads[i] && reads[i]->next && reads[i]->next->next){
				mergeSort(reads[i], readNumberByChr[i]);
			}

			size_t local_readNum = localReadNumberByChr[i];

			reads[i] = reads[i]->next;

			size_t *local_reads_coordinates_unsorted; //the vector we sort
			local_reads_coordinates_unsorted = (size_t*)malloc(local_readNum*sizeof(size_t));
			local_reads_coordinates_unsorted[0] = 0;

			size_t *local_reads_coordinates_sorted; //the vector we sort
			local_reads_coordinates_sorted = (size_t*)malloc(local_readNum*sizeof(size_t));
			local_reads_coordinates_sorted[0] = 0;

			int *local_reads_sizes_unsorted; //vector of read length
			local_reads_sizes_unsorted = (int*)malloc(local_readNum*sizeof(int));
			local_reads_sizes_unsorted[0] = 0;

			int *local_reads_sizes_sorted; //vector of read length
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
			// and indexed pbs is for parallel bitonic sort
			size_t *pbs_local_reads_coordinates;
			size_t *pbs_global_reads_coordinates_index;
			size_t dataSize = 0;

			//task Init offset and size for source - free chr
			// from mpiSort_utils.c
			get_coordinates_and_offset_source_and_size_and_free_reads(split_rank, local_reads_rank_unsorted,
					local_reads_coordinates_unsorted, local_offset_source_unsorted, local_reads_sizes_unsorted,
					reads[i], local_readNum);

			//init indices for qksort
			size_t *coord_index = (size_t*)malloc(local_readNum*sizeof(size_t));

			for(j = 0; j < local_readNum; j++){
				coord_index[j] = j;
			}

			//qksort
			base_arr2 = local_reads_coordinates_unsorted;
			qksort(coord_index, local_readNum, sizeof(size_t), 0, local_readNum - 1, compare_size_t);

			//We index data
			for(j = 0; j < local_readNum; j++){
				local_reads_coordinates_sorted[j] = local_reads_coordinates_unsorted[coord_index[j]];
				local_reads_rank_sorted[j] = local_reads_rank_unsorted[coord_index[j]];
				local_reads_sizes_sorted[j] = local_reads_sizes_unsorted[coord_index[j]];
				local_offset_source_sorted[j] = local_offset_source_unsorted[coord_index[j]];
			}



			free(coord_index); //ok
			free(local_reads_rank_unsorted); //ok
			free(local_reads_coordinates_unsorted); //ok
			free(local_reads_sizes_unsorted); //ok
			free(local_offset_source_unsorted); //ok

			MPI_Barrier(split_comm);

			//first we get the vector of all offset in destination file
			size_t *all_reads_coordinates = NULL;
			size_t *all_offsets_sources = NULL;
			int *all_reads_sizes = NULL;
			int *all_reads_rank = NULL;
			size_t *all_reads_coordinates_index = NULL;

			size_t total_num_read = 0;

			MPI_Reduce(&localReadNumberByChr[i], &total_num_read, 1, MPI_LONG_LONG_INT, MPI_SUM, chosen_split_rank, split_comm);

			if (split_rank == chosen_split_rank)
						fprintf(stderr,	"rank %d :::::[mpiSort] total_num_read = %zu \n", split_rank, total_num_read);


			MPI_Barrier(split_comm);

			if (split_rank == chosen_split_rank){

				all_reads_coordinates = (size_t *) malloc (total_num_read * sizeof(size_t));
				all_offsets_sources = (size_t *) malloc (total_num_read * sizeof(size_t));
				all_reads_sizes = (int *) malloc (total_num_read * sizeof(int));
				all_reads_rank = (int *) malloc (total_num_read * sizeof(int));

				size_t k =0;

				for (k = 0; k < total_num_read; k++){
					all_reads_sizes[k] = 0;
					all_reads_coordinates[k] = 0;
					all_offsets_sources[k] = 0;
					all_reads_rank[k] = 0;
				}
			}
			// we initialize offset_dest, the final vector
			// that is dispatch to all jobs

			MPI_Barrier(split_comm);
			// we broadcast the total number of reads to each rank
			MPI_Bcast(&total_num_read, 1, MPI_LONG_LONG_INT, chosen_split_rank, split_comm );
			//fprintf(stderr, "Rank %d ::::: after broadcast total read number = %zu \n",	split_rank, total_num_read);

			// vector of number of read per jobs
			size_t *num_reads_per_jobs = (size_t *) malloc(split_size * sizeof(size_t));
			// start_num_reads_per_jobs is the start index when dispatching the total
			size_t *start_num_reads_per_jobs = (size_t *) malloc((split_size + 1) * sizeof(size_t));

			// chosen_rank recieves the number
			// of reads of each rank and put it  in a vector
			MPI_Gather(&local_readNum, 1, MPI_LONG_LONG_INT, &num_reads_per_jobs[split_rank - chosen_split_rank], 1,
					MPI_LONG_LONG_INT, chosen_split_rank , split_comm);

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
			 * split_chosen_rank
			 * Collect sizes, coordinates and offsets
			 * in all_vector
			 */
			time_count = MPI_Wtime();
			time_count1 = MPI_Wtime();
			if (split_rank ==chosen_split_rank){

				MPI_Status status;
				//we copy the first elements in
				int k=0;
				size_t st = start_num_reads_per_jobs[chosen_split_rank];
				for (k = 0; k < num_reads_per_jobs[chosen_split_rank]; k++){
					all_reads_rank[st] = local_reads_rank_sorted[k];
					all_reads_sizes[st] = local_reads_sizes_sorted[k];
					all_offsets_sources[st] = local_offset_source_sorted[k];
					all_reads_coordinates[st] = local_reads_coordinates_sorted[k];
					st++;
				}

				for(j = 0; j < split_size; j++){
					if(j != chosen_split_rank){

						// first we care for ranks
						int *temp_buf =(int *) malloc(num_reads_per_jobs[j]* sizeof(int));
						int *temp_buf1 =(int *) malloc(num_reads_per_jobs[j]* sizeof(int));
						size_t *temp_buf2 =(size_t *) malloc(num_reads_per_jobs[j]* sizeof(size_t));
						size_t *temp_buf3 =(size_t *) malloc(num_reads_per_jobs[j]* sizeof(size_t));

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
				fprintf(stderr,	"rank %d :::::[mpiSort] Time to collect before bitonic part  %f seconds\n", split_rank, MPI_Wtime() - time_count);

			FREE_IF_NOT_NULL(local_reads_sizes_sorted); //ok
			FREE_IF_NOT_NULL(local_reads_rank_sorted); //ok
			FREE_IF_NOT_NULL(local_offset_source_sorted); //ok
			FREE_IF_NOT_NULL(local_reads_coordinates_sorted); //ok

			/*
			 * ENTER BITONIC PART
			 *
			 * The reads are sorted according
			 * to their coordinates.
			 *
			 */
			// the master rank compute the number of
			// dimension is the number of processors where we
			// perform the bitonic sort
			// int dimensions = (int)(log2(num_processes));
			// find next ( must be greater) power, and go one back
			// we compute dimension for the parallele
			// bitonic sort
			int dimensions = 1;
			while (dimensions <= split_size)
				dimensions <<= 1;

			dimensions >>= 1;

			//fprintf(stderr, "Rank %d :::::[mpiSORT] total_num_read = %zu \n", split_rank, total_num_read);
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

			MPI_Barrier(split_comm);

			if (split_rank < dimensions){

				// pbs_local_reads_coordinates is a table containing the unsorted
				// reference coordinates
				pbs_local_reads_coordinates = (size_t *)malloc(sizeof(size_t) * pbs_local_num_read_per_job[split_rank]);
				//pbs_local_reads_sizes = (int *)malloc(sizeof(int) * pbs_local_num_read_per_job[split_rank]);

				//now the master send
				//fprintf(stderr, "Rank %d :::::[WRITE] dispatch of reads coordinates for bitonic dimensions = %d \n", split_rank, dimensions);
				time_count = MPI_Wtime();
				if ( split_rank != chosen_split_rank ){

						MPI_Status status;
						//fprintf(stderr, "Rank %d :::::[WRITE] MPI_Recv from %d\n", rank, master_job_phase_1);
						MPI_Recv(pbs_local_reads_coordinates, pbs_local_num_read_per_job[split_rank],
									MPI_LONG_LONG_INT, chosen_split_rank, 0, split_comm, &status);

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

							MPI_Send(&all_reads_coordinates[pbs_start_num_coordinates_per_jobs[j]],
									pbs_local_num_read_per_job[j], MPI_LONG_LONG_INT, j, chosen_split_rank, split_comm);
						}
					}
				}
				if (split_rank == chosen_split_rank)
					fprintf(stderr,	"rank %d :::::[mpiSort] Time to dispatch all_reads_coordinates  %f seconds\n", split_rank, MPI_Wtime() - time_count);


				if (split_rank == chosen_split_rank){
					FREE_IF_NOT_NULL(all_reads_coordinates);
				}


				// we build pbs_local_reads_coordinates_index
				pbs_global_reads_coordinates_index = (size_t *)malloc(pbs_local_num_read_per_job[split_rank]*sizeof(size_t));

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
					fprintf(stderr,	"rank %d :::::[mpiSort] Call bitonic with dimensions = %d \n", split_rank, dimensions);

				ParallelBitonicSort(split_comm, split_rank, dimensions, pbs_local_reads_coordinates,
						pbs_global_reads_coordinates_index, pbs_local_num_read_per_job[split_rank],
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

				while (pbs_local_reads_coordinates[start_index] == 0){
					start_index++;
				}

				FREE_IF_NOT_NULL(pbs_local_reads_coordinates); //ok

				pbs_local_num_read_per_job[split_rank] -=  start_index;
				all_reads_coordinates_index = (size_t *)malloc(sizeof(size_t) * total_num_read);

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

				chosen_split_rank_gather_size_t(split_comm, split_rank, dimensions, chosen_split_rank,
						pbs_local_num_read_per_job[split_rank],
							pbs_local_num_read_per_job, pbs_start_num_coordinates_per_jobs,
								all_reads_coordinates_index, pbs_global_reads_coordinates_index, start_index);


				if (split_rank == chosen_split_rank)
					fprintf(stderr,	"rank %d :::::[mpiSort] Time to gather all_reads_coordinates_index %f seconds\n", split_rank, MPI_Wtime() - time_count);

				FREE_IF_NOT_NULL(pbs_global_reads_coordinates_index); //ok
				FREE_IF_NOT_NULL(pbs_local_num_read_per_job);

			} //end if (split_rank < dimensions)

			/*
			 * Bitonic is finished
			 *
			 * Now we re-order and dispatch
			 * the ranks, source offset and sizes
			 *
			 */

			FREE_IF_NOT_NULL(pbs_start_num_coordinates_per_jobs); //ok

			// we create a datatype by split_rank
			size_t *all_offset_dest_sorted=NULL;
			int *all_reads_size_sorted = NULL;
			int *all_reads_rank_sorted = NULL;
			size_t *all_offsets_sources_sorted = NULL;

			if (split_rank == chosen_split_rank){

				size_t k;
				//we reorder all_reads_sizes
				all_reads_size_sorted = (int *)malloc(sizeof(int) * total_num_read);
				all_reads_rank_sorted= (int *)malloc(sizeof(int) * total_num_read);
				all_offsets_sources_sorted = (size_t *)malloc(sizeof(size_t) * total_num_read);

				for (k = 0; k < total_num_read ; k++){
					all_reads_size_sorted[k] = all_reads_sizes[all_reads_coordinates_index[k]];
					all_offsets_sources_sorted[k] = all_offsets_sources[all_reads_coordinates_index[k]];
					all_reads_rank_sorted[k] = all_reads_rank[all_reads_coordinates_index[k]];

				}

				FREE_IF_NOT_NULL(all_reads_sizes); //ok
				FREE_IF_NOT_NULL(all_offsets_sources); //ok
				FREE_IF_NOT_NULL(all_reads_rank);

				// all_offset_dest holds the output offset of sorted reads
				all_offset_dest_sorted = (size_t *)malloc(sizeof(size_t) * total_num_read);
				all_offset_dest_sorted[0] = all_reads_size_sorted[0];
				all_offset_dest_sorted[0] += headerSize;

				for (k = 1; k < total_num_read ; k++){
					all_offset_dest_sorted[k] = all_offset_dest_sorted[k - 1] + all_reads_size_sorted[k];
				}


			} //end if (split_rank == chosen_rank)

			size_t *local_dest_offsets_sorted = (size_t*)malloc(sizeof(size_t)*local_readNum);
			size_t *local_source_offsets_sorted = (size_t*)malloc(sizeof(size_t)*local_readNum);
			int *local_read_size_sorted = (int*)malloc(sizeof(size_t)*local_readNum);
			int *local_rank_sorted = (int*)malloc(sizeof(size_t)*local_readNum);

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

					local_dest_offsets_sorted[k] = all_offset_dest_sorted[ind];
					local_source_offsets_sorted[k] = all_offsets_sources_sorted[ind];
					local_read_size_sorted[k] = all_reads_size_sorted[ind];
					local_rank_sorted[k] = all_reads_rank_sorted[ind];

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
				fprintf(stderr,	"rank %d :::::[mpiSort] Time dispatch all_offset_dest %f seconds\n", split_rank, MPI_Wtime() - time_count);

			FREE_IF_NOT_NULL(start_num_reads_per_jobs); // ok
			FREE_IF_NOT_NULL(num_reads_per_jobs); // ok


			if (split_rank == chosen_split_rank){
				FREE_IF_NOT_NULL(all_reads_rank_sorted);
				FREE_IF_NOT_NULL(all_offsets_sources_sorted);
				FREE_IF_NOT_NULL(all_offset_dest_sorted);
				FREE_IF_NOT_NULL(all_reads_size_sorted);
			}


			MPI_Barrier(split_comm);

			if (split_rank < dimensions){
				FREE_IF_NOT_NULL(all_reads_coordinates_index); //ok
			}

			writeSam(total_num_read, chosen_split_rank, split_size, dimensions, output_dir, file_name, mpi_file_split_comm, finfo, header, chrNames[i], reads[i],
					local_dest_offsets_sorted, local_source_offsets_sorted, local_read_size_sorted, local_rank_sorted,
						localReadNumberByChr[i], compression_level, split_comm);

			if (split_rank == chosen_split_rank){
				fprintf(stderr,	"rank %d :::::[mpiSort] Time spend writeSam chromosom %s ,  %f seconds\n", split_rank, chrNames[i], MPI_Wtime() - time_count);
				fprintf(stderr,	"rank %d :::::[mpiSort] Time phase1 %f seconds\n", split_rank, MPI_Wtime() - time_count1);

			}

			//FREE_IF_NOT_NULL(local_source_offsets_sorted);
			//FREE_IF_NOT_NULL(local_read_size_sorted);
			//FREE_IF_NOT_NULL(local_dest_offsets_sorted);
			//FREE_IF_NOT_NULL(local_rank_sorted);


		} //if ((local_color == 0) && (i < (nbchr - 2))) //in the splitted dimension
		else{
			//rank 0 gather the vector
			//we do nothing in this
		}

		//we put a barrier before freeing pointers
		MPI_Barrier(MPI_COMM_WORLD);
		//we free the file pointer
		if  (file_pointer_to_free)
			MPI_File_close(&mpi_file_split_comm);
		//we free the split_comm
		if (split_comm_to_free){
			MPI_Comm_free(&split_comm);
		}

		free(color_vec_to_send);
		free(key_vec_to_send);

	}// end loop upon chromosoms


	/*
	 *
	 * NOW we care for unmapped and discordants reads
	 *
	 */


	// first we test if the there's reads to sort
	// rank 0 recieve the sum of all the reads count

	int s=0;
	for (s=1; s < 3; s++){

		MPI_File mpi_file_split_comm2;
		double time_count;

		size_t total_reads = 0;
		MPI_Allreduce(&readNumberByChr[nbchr-s], &total_reads , 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
		if ((rank == 0) && (s == 1))
			fprintf(stderr, "rank %d :::: total read to sort for unmapped = %zu \n", rank, total_reads);

		if ((rank == 0) && (s == 2))
			fprintf(stderr, "rank %d :::: total read to sort for discordant = %zu \n", rank, total_reads);

		MPI_Barrier(MPI_COMM_WORLD);

		if (total_reads == 0){
			// nothing to sort for unmapped
			// maybe write an empty bam file
			// we go directly to discordant reads
		}
		else{

			int i1,i2;
			size_t *localReadsNum_rank0 = (size_t *)malloc(num_proc*sizeof(size_t));
			int file_pointer_to_free = 0;
			int split_comm_to_free = 0;
			//we build a vector with rank job
			int val_tmp1 = 0;
			int val_tmp2 = 0;
			int chosen_rank = 0;
			// the color tells in what communicator the rank pertain
			// color = 0 will be the new communicator color
			// otherwise the color is 1
			int *color_vec_to_send =  (int *)malloc(num_proc*sizeof(int));
			// the key value tell the order in the new communicator
			int *key_vec_to_send =  (int *)malloc(num_proc*sizeof(int));

			//rank 0 gather the vector
			MPI_Allgather(&readNumberByChr[nbchr-s] , 1, MPI_LONG_LONG_INT, localReadsNum_rank0 , 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);

			if (rank == 0){
				//we must chose the first rank with reads to sort
				i1=0;
				while (localReadsNum_rank0[i1] == 0){
					chosen_rank++;
					i1++;
				}
			}

			//we broadcast the chosen rank
			//task: replace the broadcast with a sendrecieve
			MPI_Bcast( &chosen_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);

			//we must chose which rank is going to split the communication
			if (((rank == chosen_rank) || rank == 0) && (chosen_rank != 0)){
				//the rank 0 will recieve the key_vec_to_send and colorvec_to_send


				//first we exchange the size o
				if (rank == chosen_rank){
					header=(char *)malloc((headerSize + 1)*sizeof(char));
					MPI_Recv(header, headerSize + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				if (rank == 0){
					MPI_Send(header, headerSize + 1, MPI_CHAR, chosen_rank,  0, MPI_COMM_WORLD);
				}
			}
			else {
				//we do nothing here
			}

			if (rank == chosen_rank) {

				int counter = 0;
				//we compute the number of 0 in the localReadsNum_vec
				for(i1 = 0; i1 < num_proc; i1++){
					if (localReadsNum_rank0[i1] == 0) {
						counter++;
					}
				}
				// if no jobs without reads we do nothing
				if ( counter == 0 ){
					// nothing to do we associate split_comm with
					split_comm = MPI_COMM_WORLD;
					for (i2 = 0; i2 < num_proc; i2++) {

						if (localReadsNum_rank0[i] == 0) {
							color_vec_to_send[i2] = 1;
							key_vec_to_send[i2] = val_tmp2;
							val_tmp2++;
						} else {
							color_vec_to_send[i2] = 0;
							key_vec_to_send[i2] = val_tmp1;
							val_tmp1++;
						}
					}
				}
				else{

					// now we compute the color according to
					// the number of reads to sort
					for(i2 = 0; i2 < num_proc; i2++){

						if (localReadsNum_rank0[i2] == 0){
							color_vec_to_send[i2] = 1;
							key_vec_to_send[i2] = val_tmp2;
							val_tmp2++;
						} else{
							color_vec_to_send[i2] = 0;
							key_vec_to_send[i2] = val_tmp1;
							val_tmp1++;
						}
					} // end for loop
				}// end if

			}// end if (rank == chosen_rank)

			MPI_Barrier(MPI_COMM_WORLD);
			// we scatter the key and color vector
			//we create key and color variable for each job
			int local_color = 0;
			int local_key = 0;
			// we scatter the color and key
			MPI_Scatter( color_vec_to_send, 1, MPI_INT, &local_color, 1, MPI_INT, chosen_rank, MPI_COMM_WORLD);
			MPI_Scatter( key_vec_to_send, 1, MPI_INT, &local_key, 1, MPI_INT, chosen_rank, MPI_COMM_WORLD);
			// now we create a communicator
			// we group all communicator
			// with color of zero

			if (local_color == 0){

				MPI_Comm_split( MPI_COMM_WORLD, local_color, local_key, &split_comm);
				ierr = MPI_File_open(split_comm, file_name,  MPI_MODE_RDONLY , finfo, &mpi_file_split_comm2);
				//we ask to liberate file pointer
				file_pointer_to_free = 1;
				//we ask to liberate the split_comm
				split_comm_to_free = 1;
			}
			else{
				MPI_Comm_split( MPI_COMM_WORLD, MPI_UNDEFINED, local_key, &split_comm);
				mpi_file_split_comm2 = mpi_filed;
			}

			//we broadcast to_free at the end of the loop
			//MPI_Bcast(&to_free, 1, MPI_INT, 0, MPI_COMM_WORLD);

			//now we change the rank in the reads structure
			if (local_color == 0){
				MPI_Comm_rank(split_comm, &split_rank);
				MPI_Comm_size(split_comm, &split_size);

				g_rank = split_rank;
				g_size = split_size;

				reads[nbchr-s] = reads[nbchr-s]->next;
				localReadNumberByChr[nbchr-s] = readNumberByChr[nbchr-s];
				if (s == 2){
					unmapped_start = startOffset(g_rank, g_size,
										unmappedSize, headerSize,
											nbchr-s, localReadNumberByChr[nbchr-s], split_comm);

					if(!unmapped_start){
						fprintf(stderr, "No header was defined for unmapped. \n Shutting down.\n");
						MPI_Finalize();
						return 0;
					}

					time_count = MPI_Wtime();

					writeSam_discordant_and_unmapped(split_rank, output_dir, header, localReadNumberByChr[nbchr-s], chrNames[nbchr-s], reads[nbchr-s],
													split_size, split_comm, file_name, mpi_file_split_comm2, finfo, compression_level);
					if (split_rank == chosen_rank){
							fprintf(stderr,	"rank %d :::::[mpiSort] Time to write chromosom %s ,  %f seconds\n", split_rank,
									chrNames[nbchr-s], MPI_Wtime() - time_count);
					}

					/*
					 *
					 * TODO in future release replace writeSam_unmapped with: writeSam_notSorted
					 *
					 * The troubles are in the global variables
					 *
					 *
						writeSam_notSorted(output_dir, header,
							localReadNumberByChr[nbchr-s], chrNames[nbchr-s],
								reads[nbchr-s], file_name,
									mpi_file_split_comm2, finfo, compression_level);
					*/



				}
				else{
					discordant_start = startOffset(g_rank, g_size,
													discordantSize, headerSize,
														nbchr-s, localReadNumberByChr[nbchr-s], split_comm);
					if(!discordant_start){
						fprintf(stderr, "No header was defined for discordant.\n Shutting down.\n");
						MPI_Finalize();
						return 0;
					}

					time_count = MPI_Wtime();

					writeSam_discordant_and_unmapped(g_rank, output_dir, header, localReadNumberByChr[nbchr-s], chrNames[nbchr-s], reads[nbchr-s],
							g_size, split_comm, file_name, mpi_file_split_comm2, finfo, compression_level);

					if (split_rank == chosen_rank){
							fprintf(stderr,	"rank %d :::::[mpiSort] Time to write chromosom %s ,  %f seconds\n", split_rank,
								chrNames[nbchr-s], MPI_Wtime() - time_count);
					}

					/*
					 *
					 * TODO in future release replace writeSam_unmapped with: writeSam_notSorted
					 *
					 * The troubles are in the global variables
					 *
					 *
					writeSam_notSorted(output_dir, header,
							localReadNumberByChr[nbchr-s], chrNames[nbchr-s],
								reads[nbchr-s], file_name,
									mpi_file_split_comm2, finfo, compression_level);
					*/


				}
				while( reads[nbchr-s]->next != NULL){
						Read *tmp_chr = reads[nbchr-s];
						reads[nbchr-s] = reads[nbchr-s]->next;
						//FREE_IF_NOT_NULL(tmp_chr->next);
				}
				FREE_IF_NOT_NULL(localReadsNum_rank0);
			}
			else{
				// we do nothing
			}

			//we put a barrier before freeing pointers
			MPI_Barrier(MPI_COMM_WORLD);
			//we free the file pointer

			/*
			 * TODO problem freeing this pointer
			 */
			if  (file_pointer_to_free)
				MPI_File_close(&mpi_file_split_comm2);

			//we free the split_comm
			if (split_comm_to_free)
				MPI_Comm_free(&split_comm);


			split_comm_to_free = 0;
			file_pointer_to_free = 0;

			FREE_IF_NOT_NULL(color_vec_to_send);
			FREE_IF_NOT_NULL(key_vec_to_send);

		}
	} //end for (s=1; s < 3; s++){
	MPI_Barrier(MPI_COMM_WORLD);

	FREE_IF_NOT_NULL(header); //ok
	FREE_IF_NOT_NULL(localReadNumberByChr); //ok

	for(i = 0; i < nbchr; i++){
		FREE_IF_NOT_NULL(chrNames[i]);
	}


	FREE_IF_NOT_NULL(reads); //ok
	FREE_IF_NOT_NULL(readNumberByChr); //ok

	res = MPI_Finalize();
	assert(res == MPI_SUCCESS);

	return 0;
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

