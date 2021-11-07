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
     mpiSORT.c

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
#include <stdlib.h>
#include <mpi.h>

#include "compat.h"
#include "malloc.h"
#include "mergeSort.h"
#include "parser.h"
#include "preWrite.h"
#include "write.h"
#include "sortAnyDim.h"
#include "mpiSortUtils.h"
#include "writeUtils.h"
#include "qkSort.h"
#include "parallelBitonicSort2.h"
#include "parallelBitonicSort3.h"


/*
 * Constant for Lustre striping
 * To adapt depending on the file server
 * STRIPING FACTOR is the numer of servers
 * where your files are striped
 *
 * STRIPING UNIT is the size of the stripes
 *
 * no need to change the values if you don't use
 * Lustre.
 *
 * Those parameters are harmless for other FS
 *
 */
#define STRIPING_FACTOR "12"
#define STRIPING_UNIT "4194304"   // 4 MB
//#define STRIPING_UNIT "268435456"   // 256 MB
//#define STRIPING_UNIT "536870912"   // 500 MB
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

// BEGIN> FINE TUNING FINFO FOR WRITING OPERATIONS
#define NB_PROC  "20" //numer of threads for writing
#define CB_NODES "12" //numer of server for writing
#define CB_BLOCK_SIZE  "268435456" /* 256 MBytes - should match FS block size */
#define CB_BUFFER_SIZE  "536870912" /* multiple of the block size by the number of proc*/
#define DATA_SIEVING_READ "enable"
// END> FINE TUNING FINFO FOR WRITING OPERATIONS

/*
 * Capacity constant no need to change it
 */
#define DEFAULT_MAX_SIZE 6000000000 //Default capacity for one process: 6G
#define DEFAULT_INBUF_SIZE  (1024*1024*1024)

/* Maximum chromosome number */
#define MAXNBCHR 512

static void usage(const char *);

int main (int argc, char *argv[]){

	char *x, *x1, *y, *y1, *z, *z1, *xbuf, *hbuf, *chrNames[MAXNBCHR];
	int fd;
	off_t hsiz;
	struct stat st;

	MPI_File mpi_filed;
	MPI_File mpi_file_split_comm;

	MPI_Offset fileSize, unmapped_start, discordant_start;
	int num_proc, rank;
	int res, nbchr, i, paired, uniq_chr, write_format;
	int ierr, errorcode = MPI_ERR_OTHER;
	char *file_name, *output_dir;

	char *header;

	unsigned int headerSize;
	unsigned char threshold;

	size_t input_file_size;
	size_t unmappedSize = 0;
	size_t discordantSize = 0;
	size_t *readNumberByChr = NULL, *localReadNumberByChr = NULL;
	Read **reads;

	double time_count;
	double time_count1;
	int g_rank, g_size;
	MPI_Comm split_comm; //used to split communication when jobs have no reads to sort
	int split_rank, split_size; //after split communication we update the rank and the size
	double tic, toc;
	int compression_level;
	int merge;
	size_t fsiz, lsiz, loff;
	const char *sort_name;
        char *chr_name_u; //use is case of uniq chromosom. Must be the name of the chromosom, it gives also the name of output file 
	MPI_Info finfo;

	/* Set default values */
	compression_level = 3;
	parse_mode = MODE_OFFSET;
	sort_name = "coordinate";
	paired = 0; /* by default reads are considered single*/
	uniq_chr = 0; 
	threshold = 0;
	write_format = 0;
	merge = 0;
	/* Check command line */
	while ((i = getopt(argc, argv, "c:hnpu:q:gsbm")) != -1) {
		switch(i) {
			case 'c': /* Compression level */
				compression_level = atoi(optarg);
				break;
			case 'h': /* Usage display */
				usage(basename(*argv));
				return 0;
			case 'n':
				parse_mode = MODE_NAME;
				sort_name = "queryname";
				break;
			case 'p': /* Paired reads */
				paired = 1;
				break;
			case 'u': /* We say we have only one chromosome in the file */
                                uniq_chr = 1;
				asprintf(&chr_name_u,"%s", optarg);
                                break;
			case 'q': /* Quality threshold */
				threshold = atoi(optarg);
				break;
			case 'g':
				write_format = 0;
				break;
			case 'b':
                                write_format = 1;
                                break;
			case 's':
                                write_format = 2;
                                break;
			case 'm':
				merge = 1;
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

	char *file_name_tmp;
        char file_name_merge[256];
        char *dot;
	char *slash = strdup(file_name);
        //in case of merge we create a file from file_name
        if (merge){

	   	file_name_tmp = basename(slash);
                dot = strrchr(file_name_tmp, '.');
                if (dot) dot[0] = '\0';
                sprintf(file_name_merge, "%s_sorted.sam", file_name_tmp);
        }

	/* MPI inits */
	res = MPI_Init(&argc, &argv);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	assert(res == MPI_SUCCESS);


	/* Check if num_proc is a power of 2 or not null */
        if (num_proc == 0){
		fprintf(stderr, "Number of processes must be greater than 0 \n");
		res = MPI_Finalize();
                assert(res == MPI_SUCCESS);
                exit(2);
                //err(1, "You ask for 0 cpu this is not possible !!\n");	
	}

	
	//remove this condition to play with non power of 2
	
        if ( (num_proc & (num_proc - 1)) ){
		fprintf(stderr, "Number of processes must be power of two \n");
                res = MPI_Finalize();
                assert(res == MPI_SUCCESS);
                exit(2);
                //err(1, "Number of processes must be power of two \n");
	}

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
	xbuf = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_FILE|MAP_PRIVATE, fd, 0);
	assert(xbuf != MAP_FAILED);
	/* we get the size of the original header */
	x = xbuf;
	size_t orig_hsiz = 0;
	while (*x == '@') {
                y = strchr(x, '\n');
                z = x; x = y + 1;
              
        }
	orig_hsiz = x - xbuf;

	/* ignore the first line if equal to @HD     VN:1.0  SO:  */
        char *xbuf1;
        x1 = xbuf;
	xbuf1 = xbuf;
	char *y2;
	
	if (*x1 == '@'){
		y1 = strchr(x1, '\n');		
		z1 = x1; x1 = y1 + 1;
		if (strncmp(z1, "@HD", 3) == 0) {
			y2 = strstr(z1, "SO:");
			if (y2 != NULL) xbuf1 = x1;	
		}					
	}
	
	/* Parse SAM header */
	memset(chrNames, 0, sizeof(chrNames));
	x = xbuf1; nbchr = 0;
	
	char* ts1 = strdup(file_name);
	char* ts2 = strdup(file_name);

        char* dir = dirname(ts1);
        char* filename = basename(ts2);
	
	if (uniq_chr){

		//char *tmp_str = filename;
		//const char *dot = strrchr(filename, '.');
		//chrNames[nbchr++] = strndup(tmp_str, dot-tmp_str);		
		chrNames[nbchr++] = strdup(chr_name_u);
	}
	
	assert(*x == '@');
	while (*x == '@') {
		y = strchr(x, '\n');
		z = x; x = y + 1;
		if (strncmp(z, "@SQ", 3) != 0) continue;
		/* Save reference names */
		y = strstr(z, "SN:");
		assert(y != NULL);
		z = y + 3;
		while (*z && !isspace((unsigned char)*z)) z++;
		if (!uniq_chr) chrNames[nbchr++] = strndup(y + 3, z - y - 3);
		assert(nbchr < MAXNBCHR - 2);
	}	
	assert(*x != '@');	
	
	
	
	//in the case of a unique chromosome in the sam
	//the discordant file is named chrX_discordant
	if (uniq_chr) {
		asprintf(&chrNames[nbchr++],"%s_%s", chr_name_u, DISCORDANT);
		asprintf(&chrNames[nbchr++],"%s_%s", chr_name_u, UNMAPPED);
	}
	else {
		chrNames[nbchr++] = strdup(DISCORDANT);
		chrNames[nbchr++] = strdup(UNMAPPED);
	}

	hsiz = x - xbuf1;
	hbuf = strndup(xbuf1, hsiz);

	if (rank == 0) {
		fprintf(stderr, "The size of the file is %zu bytes\n", (size_t)st.st_size);
		fprintf(stderr, "Header has %d+2 references\n", nbchr - 2);
	}
	asprintf(&header, "@HD\tVN:1.0\tSO:%s\n%s", sort_name, hbuf);

	free(hbuf);

	assert(munmap(xbuf, (size_t)st.st_size) != -1);
	assert(close(fd) != -1);


    // BEGIN> FINE TUNING FINFO FOR WRITING OPERATIONS

	MPI_Info_create(&finfo);
	/*
	 * In this part you shall adjust the striping factor and unit according
	 * to the underlying filesystem.
	 * Harmless for other file system.
	 *
	 */
	MPI_Info_set(finfo,"striping_factor", STRIPING_FACTOR);
	MPI_Info_set(finfo,"striping_unit", STRIPING_UNIT); //2G striping
	MPI_Info_set(finfo,"ind_rd_buffer_size", STRIPING_UNIT); //2gb buffer
	MPI_Info_set(finfo,"romio_ds_read",DATA_SIEVING_READ);

	/*
	 * for collective reading and writing
	 * should be adapted too and tested according to the file system
	 * Harmless for other file system.
	 */
	MPI_Info_set(finfo,"nb_proc", NB_PROC);
	MPI_Info_set(finfo,"cb_nodes", CB_NODES);
	MPI_Info_set(finfo,"cb_block_size", CB_BLOCK_SIZE);
	MPI_Info_set(finfo,"cb_buffer_size", CB_BUFFER_SIZE);
    
    // END> FINE TUNING FINFO FOR WRITING OPERATIONS


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

	/* Get chunk offset and size */
	fsiz = input_file_size;
	lsiz = fsiz / num_proc;
	loff = rank * lsiz;

	tic = MPI_Wtime();

	headerSize = unmappedSize = discordantSize = strlen(header);

	//We place file offset of each process to the begining of one read's line
	size_t *goff =(size_t*)calloc((size_t)(num_proc+1), sizeof(size_t));
	init_goff(mpi_filed,orig_hsiz,input_file_size,num_proc,rank,goff);

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

	char *local_data_tmp = malloc(1024*1024);
	char *local_data =(char*)malloc(((goff[rank+1]-poffset)+1)*sizeof(char));
	size_t size_tmp= goff[rank+1]-poffset;
	local_data[goff[rank+1]-poffset] = 0;
	char *q=local_data;

	//We read the file sam and parse
	while(poffset < goff[rank+1]){

		size_t size_to_read = 0;

		if( (goff[rank+1]-poffset) < DEFAULT_INBUF_SIZE ){
			size_to_read = goff[rank+1]-poffset;
		}
		else{
			size_to_read = DEFAULT_INBUF_SIZE;
		}

		// we load the buffer
		//hold temporary size of SAM
		//due to limitation in MPI_File_read_at
		local_data_tmp =(char*)realloc(local_data_tmp, (size_to_read+1)*sizeof(char));
		local_data_tmp[size_to_read]=0;

		// Original reading part is before 18/09/2015
		MPI_File_read_at(mpi_filed, (MPI_Offset)poffset, local_data_tmp, size_to_read, MPI_CHAR, MPI_STATUS_IGNORE);
		size_t local_offset=0;
		assert(strlen(local_data_tmp) == size_to_read);

		//we look where is the last line read for updating next poffset
		size_t offset_last_line = size_to_read-1;

		size_t extra_char=0;
		while(local_data_tmp[offset_last_line] != '\n'){
			offset_last_line -- ;
			extra_char++;
		}

		local_data_tmp[size_to_read - extra_char]=0;
		size_t local_data_tmp_sz = strlen(local_data_tmp);

		//If it s the last line of file, we place a last '\n' for the function tokenizer
		if(rank == num_proc-1 && ((poffset+size_to_read) == goff[num_proc])){
			local_data_tmp[offset_last_line]='\n';
		}

		//Now we parse Read in local_data
		if (paired == 1 && uniq_chr == 0) parser_paired(local_data_tmp, rank, poffset, threshold, nbchr, &readNumberByChr, chrNames, &reads);
		if (paired == 1 && uniq_chr == 1) parser_paired_uniq(local_data_tmp, rank, poffset, threshold, nbchr, &readNumberByChr, chrNames, &reads);
		if (paired == 0) parser_single(local_data_tmp, rank, poffset, threshold, nbchr, &readNumberByChr, chrNames, &reads);
		//now we copy local_data_tmp in local_data
		char *p = local_data_tmp;
		int pos =0;
		while (*p && (pos < local_data_tmp_sz)) {*q=*p;p++;q++;pos++;}

		//we go to the next line
		poffset+=(offset_last_line+1);
		local_offset+=(offset_last_line+1);

	}

	assert(size_tmp == strlen(local_data));

	fprintf(stderr, "%d (%.2lf)::::: *** FINISH PARSING FILE ***\n", rank, MPI_Wtime()-toc);

	if (local_data_tmp) free(local_data_tmp);
	malloc_trim(0);

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
	if (rank == 0)
	 	fprintf(stderr, "rank %d ::::[MPISORT] total reads parsed = %zu \n", rank, nb_reads_global);
	/*
	 * We care for unmapped and discordants reads
	 */

	int s = 0;
	for (s = 1; s < 3; s++){

		MPI_File mpi_file_split_comm2;
		double time_count;

		size_t total_reads = 0;
		MPI_Allreduce(&readNumberByChr[nbchr-s], &total_reads , 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

		if ((rank == 0) && (s == 1))
			fprintf(stderr, "rank %d ::::[MPISORT] total read to sort for unmapped = %zu \n", rank, total_reads);

		if ((rank == 0) && (s == 2))
			fprintf(stderr, "rank %d ::::[MPISORT] total read to sort for discordant = %zu \n", rank, total_reads);

		MPI_Barrier(MPI_COMM_WORLD);

		if (total_reads == 0){
			// nothing to sort for unmapped
			// maybe write an empty bam file
		}
		else{
			int i1,i2;
			size_t *localReadsNum_rank0 = (size_t *)malloc(num_proc*sizeof(size_t));
			localReadsNum_rank0[0] = 0;
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

						if (localReadsNum_rank0[i2] == 0) {
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
			// we create key and color variable for each job
			int local_color = 0;
			int local_key = 0;
			// we scatter the color and key
			MPI_Scatter( color_vec_to_send, 1, MPI_INT, &local_color, 1, MPI_INT, chosen_rank, MPI_COMM_WORLD);
			MPI_Scatter( key_vec_to_send, 1, MPI_INT, &local_key, 1, MPI_INT, chosen_rank, MPI_COMM_WORLD);
			// we create a communicator
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

			//now we change the rank in the reads structure
			if (local_color == 0){
				MPI_Comm_rank(split_comm, &split_rank);
				MPI_Comm_size(split_comm, &split_size);

				g_rank = split_rank;
				g_size = split_size;

				reads[nbchr-s] = reads[nbchr-s]->next;
				localReadNumberByChr[nbchr-s] = readNumberByChr[nbchr-s];
				if (s == 2){
					unmapped_start = startOffset(g_rank,
												 g_size,
												 unmappedSize,
												 headerSize,
												 nbchr-s,
												 localReadNumberByChr[nbchr-s],
												 split_comm
												 );

					if(!unmapped_start){
						fprintf(stderr, "No header was defined for unmapped. \n Shutting down.\n");
						MPI_Finalize();
						return 0;
					}

					time_count = MPI_Wtime();
					writeSam_discordant_and_unmapped(
							split_rank,
							output_dir,
							header,
							localReadNumberByChr[nbchr-s],
							chrNames[nbchr-s],
							reads[nbchr-s],
							split_size,
							split_comm,
							file_name,
							mpi_file_split_comm2,
							finfo,
							compression_level,
							local_data,
							goff[rank],
							write_format);

					if (split_rank == chosen_rank){
							fprintf(stderr,	"rank %d :::::[MPISORT] Time to write chromosom %s ,  %f seconds \n\n\n", split_rank,
									chrNames[nbchr-s], MPI_Wtime() - time_count);
					}
				}
				else{
					discordant_start = startOffset(g_rank,
												   g_size,
												   discordantSize,
												   headerSize,
												   nbchr-s,
												   localReadNumberByChr[nbchr-s],
												   split_comm);

					if(!discordant_start){
						fprintf(stderr, "No header was defined for discordant.\n Shutting down.\n");
						MPI_Finalize();
						return 0;
					}
					time_count = MPI_Wtime();

					writeSam_discordant_and_unmapped(
							g_rank,
							output_dir,
							header,
							localReadNumberByChr[nbchr-s],
							chrNames[nbchr-s],
							reads[nbchr-s],
							g_size,
							split_comm,
							file_name,
							mpi_file_split_comm2,
							finfo,
							compression_level,
							local_data,
							goff[rank],
							write_format
							);


					if (split_rank == chosen_rank){
							fprintf(stderr,	"rank %d :::::[MPISORT] Time to write chromosom %s ,  %f seconds \n\n\n", split_rank,
								chrNames[nbchr-s], MPI_Wtime() - time_count);
					}

				}
				while( reads[nbchr-s]->next != NULL){
						Read *tmp_chr = reads[nbchr-s];
						reads[nbchr-s] = reads[nbchr-s]->next;
						free(tmp_chr);
				}
				free(localReadsNum_rank0);
			}
			else{
				// we do nothing
			}

			//we put a barrier before freeing pointers
			MPI_Barrier(MPI_COMM_WORLD);
			//we free the file pointer

			if  (file_pointer_to_free)
				MPI_File_close(&mpi_file_split_comm2);

			//we free the split_comm
			if (split_comm_to_free)
				MPI_Comm_free(&split_comm);

			split_comm_to_free = 0;
			file_pointer_to_free = 0;

			free(color_vec_to_send);
			free(key_vec_to_send);

		}
	} //end for (s=1; s < 3; s++){

	/*
	 *  We write the mapped reads in a file named chrX.bam
	 *	We loop by chromosoms.
	 */

	MPI_Barrier(MPI_COMM_WORLD);
	for(i = 0; i < (nbchr-2); i++){
		/*
		 * First Part of the algorithm
		 *
		 * In this part we elected a rank which is the first rank
		 * to have reads to sort.
		 *
		 * Once elected a rank, we plit the communicator according to
		 * wether the rank has reads to sort for this chromosom.
		 *
		 * The new communicator is COMM_WORLD.
		 *
		 * If all jobs have reads to sort no need to split the communicator and then
		 * COMM_WORLD = MPI_COMM_WORLD
		 *
		 */

		int i1,i2;
		size_t localReadsNum_rank0[num_proc];
		localReadsNum_rank0[0]=0;
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
		// the key value tell the order in the new communicator
		int *color_vec_to_send 	=  malloc(num_proc * sizeof(int));
		int *key_vec_to_send 	=  malloc(num_proc * sizeof(int));

		// first we test if the there's reads to sort
		// rank 0 recieve the sum of all the reads count
		size_t total_reads_by_chr = 0;
		MPI_Allreduce(&readNumberByChr[i], &total_reads_by_chr, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

		//fprintf(stderr, "rank %d :::: readNumberByChr[i] = %zu \n", rank, readNumberByChr[i]);
		//fprintf(stderr, "rank %d :::: total_reads_by_chr = %zu \n", rank, total_reads_by_chr);

		if (total_reads_by_chr == 0)
			continue; //pass to next chromosome

		//rank 0 gather the vector
		MPI_Allgather(&readNumberByChr[i] , 1, MPI_LONG_LONG_INT, localReadsNum_rank0 , 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);


		if (rank == 0){
			//the rank 0 chose the first rank with reads to sort
			i1=0;
			while ((localReadsNum_rank0[i1] == 0) && (i1 < num_proc)){
				chosen_rank++;
				i1++;
			}
			fprintf(stderr, "rank %d ::::[MPISORT] Elected rank = %d \n", rank, chosen_rank);
		}

		//we broadcast the chosen rank
		//task: replace the broadcast with a sendrecieve
		MPI_Bcast( &chosen_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		if (((rank == chosen_rank) || rank == 0) && (chosen_rank != 0)){

			//first we exchange the size o
			if (rank == chosen_rank){
				header = malloc((headerSize + 1)*sizeof(char));
				header[headerSize] = '\0';
				MPI_Recv(header, headerSize + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			if (rank == 0){
				MPI_Send(header, headerSize + 1, MPI_CHAR, chosen_rank,  0, MPI_COMM_WORLD);
			}
		}
		else {
			//we do nothing here
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
				fprintf(stderr, "rank %d ::::[MPISORT] we don't split the rank \n", rank);
				split_comm = MPI_COMM_WORLD;
				for (i2 = 0; i2 < num_proc; i2++) {
					if (localReadsNum_rank0[i2] == 0) {
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
				fprintf(stderr, "rank %d ::::[MPISORT] we split the rank \n", rank);
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
			g_rank = split_rank;			
			g_size = split_size = num_proc;
		}

		localReadNumberByChr[i] = readNumberByChr[i];
		MPI_Barrier(MPI_COMM_WORLD);

		if ((local_color == 0) && (i < (nbchr - 2))) {

			/*
			 * Second part of the algorithm
			 *
			 * First we load coordinates, offset sources, and read size in vector
			 *
			 * Then we sort the coordinates of the reads
			 * with a bitonic sorter
			 *
			 * Then according to the reads coordinates we reoder the offset sources, and size
			 * this is done thanks to the index of the sorting.
			 *
			 * Afterward we compute the offsets of the reads in
			 * the destination file.
			 *
			 * Finally we dispatch the information to all ranks
			 * in the communicator for the next step.
			 */

			//we do a local merge sort
			if(reads[i] && reads[i]->next && reads[i]->next->next){
				mergeSort(reads[i], readNumberByChr[i]);
			}

			size_t local_readNum = localReadNumberByChr[i];

			reads[i] = reads[i]->next;

			//first we compute the dimension of the parabitonic sort
			// dimension is the number of processors where we
			// perform the bitonic sort
			// int dimensions = (int)(log2(num_processes));
			// find next ( must be greater) power, and go one back
			int dimensions = 1;
			while (dimensions <= split_size)
				dimensions <<= 1;

			dimensions >>= 1;

			// we get the maximum number of reads among
			// all the workers

			/*
			 * Here we split the programm in 2 cases
			 *
			 * 1) The first case de split_size is a power of 2 (the best case)
			 * 		this case is the simpliest we don't have extra communication to dispatch the read
			 * 		envenly between the jobs
			 *
			 * 2) The split_size is not a power of 2 (the worst case)
			 * 		well in this case we shall dispatch the jobs between jobs evenly.
			 *
			 */

			if (split_rank == chosen_split_rank){

				fprintf(stderr,	"Rank %d :::::[MPISORT] Dimensions for bitonic = %d \n", split_rank, dimensions);
				fprintf(stderr,	"Rank %d :::::[MPISORT] Split size 			   = %d \n", split_rank, split_size);

			}
			//we test the computed dimension
			if (dimensions == split_size ){

				size_t max_num_read = 0;
				MPI_Allreduce(&localReadNumberByChr[i], &max_num_read, 1, MPI_LONG_LONG_INT, MPI_MAX, split_comm);

				// if the dimension == split_size
				MPI_Barrier(split_comm);

				size_t first_local_readNum = local_readNum;

				/*
				 * Vector creation and allocation
				 				fprintf(stderr,	"split rank %d :::::[MPISORT] max_num_read = %zu \n", split_rank, max_num_read);
				 */
				local_readNum = max_num_read;

				time_count = MPI_Wtime();

				size_t *local_reads_coordinates_unsorted 	= calloc(local_readNum, sizeof(size_t));
				size_t *local_reads_coordinates_sorted 		= calloc(local_readNum, sizeof(size_t));
				size_t *local_offset_source_unsorted 		= calloc(local_readNum, sizeof(size_t));
				size_t *local_offset_source_sorted 		= calloc(local_readNum, sizeof(size_t));
				int *local_dest_rank_sorted 			= calloc(local_readNum, sizeof(int));
				int *local_reads_sizes_unsorted 		= calloc(local_readNum, sizeof(int));
				int *local_reads_sizes_sorted 			= calloc(local_readNum, sizeof(int));
				int *local_source_rank_unsorted 		= calloc(local_readNum, sizeof(int));
				int *local_source_rank_sorted 			= calloc(local_readNum, sizeof(int));

				if (split_rank == chosen_split_rank)
					fprintf(stderr,	"rank %d :::::[MPISORT][MALLOC 1] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);

				local_reads_coordinates_unsorted[0] = 0;
				local_reads_coordinates_sorted[0] 	= 0;
				local_dest_rank_sorted[0] 			= 0;
				local_reads_sizes_unsorted[0] 		= 0;
				local_reads_sizes_sorted[0] 		= 0;
				local_source_rank_unsorted[0] 		= 0;
				local_source_rank_sorted[0] 		= 0;
				local_offset_source_unsorted[0] 	= 0;
				local_offset_source_sorted[0] 		= 0;

				//those vectors are the same that  local_..._sorted but without zero padding
				size_t *local_reads_coordinates_sorted_trimmed = NULL;
				int *local_dest_rank_sorted_trimmed = NULL;
				int *local_reads_sizes_sorted_trimmed = NULL;
				size_t *local_offset_source_sorted_trimmed = NULL;
				size_t *local_offset_dest_sorted_trimmed = NULL;
				int *local_source_rank_sorted_trimmed = NULL;

				//vectors used in the bruck just after the parabitonic sort
				size_t *local_reads_coordinates_sorted_trimmed_for_bruck = NULL;
				int *local_dest_rank_sorted_trimmed_for_bruck = NULL;
				int *local_reads_sizes_sorted_trimmed_for_bruck = NULL;
				size_t *local_offset_source_sorted_trimmed_for_bruck = NULL;
				size_t *local_offset_dest_sorted_trimmed_for_bruck = NULL;
				int *local_source_rank_sorted_trimmed_for_bruck = NULL;


				//task Init offset and size for source - free chr
				// from mpiSort_utils.c
				get_coordinates_and_offset_source_and_size_and_free_reads(
						split_rank,
						local_source_rank_unsorted,
						local_reads_coordinates_unsorted,
						local_offset_source_unsorted,
						local_reads_sizes_unsorted,
						reads[i],
						first_local_readNum
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
				time_count = MPI_Wtime();

				base_arr2 = local_reads_coordinates_unsorted;
				
				// For stabilization we replace the qksort with a merge sort
				//qksort(coord_index, local_readNum, sizeof(size_t), 0, local_readNum - 1, compare_size_t);
				MergeSortMain(coord_index, local_readNum);

				if (split_rank == chosen_split_rank)
						fprintf(stderr,	"rank %d :::::[MPISORT][LOCAL SORT] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);

				//We index data
				for(j = 0; j < local_readNum; j++){
					local_reads_coordinates_sorted[j] 			= local_reads_coordinates_unsorted[coord_index[j]];
					local_source_rank_sorted[j] 				= local_source_rank_unsorted[coord_index[j]];
					local_reads_sizes_sorted[j] 				= local_reads_sizes_unsorted[coord_index[j]];
					local_offset_source_sorted[j] 				= local_offset_source_unsorted[coord_index[j]];
					local_dest_rank_sorted[j] 					= rank; //will be updated after sorting the coordinates
				}

				/*
				*   FOR DEBUG
				*/  
					

				for(j = 0; j < local_readNum - 1; j++){
					assert( local_reads_coordinates_sorted[j] <= local_reads_coordinates_sorted[j+1]);
				}
				

				free(coord_index); 				 		//ok
				free(local_source_rank_unsorted); 	    //ok
				free(local_reads_coordinates_unsorted); //ok
				free(local_reads_sizes_unsorted); 		//ok
				free(local_offset_source_unsorted); 	//ok

				// we need the total number of reads.
				size_t total_num_read = 0;
				MPI_Allreduce(&localReadNumberByChr[i], &total_num_read, 1, MPI_LONG_LONG_INT, MPI_SUM, split_comm);

				/*
				 *
				 * In this section the number of bitonic dimension
				 * is equal to the split size.
				 *
				 * In this case there are less communication in preparation
				 * of the sorting.
				 *
				 * We use the parabitonic version 2.
				 */

				//we calll the bitonic

				time_count = MPI_Wtime();

				ParallelBitonicSort2(
					split_comm,
					split_rank,
					dimensions,
					local_reads_coordinates_sorted,
					local_reads_sizes_sorted,
					local_source_rank_sorted,
					local_offset_source_sorted,
					local_dest_rank_sorted,
					max_num_read
					);

				if (split_rank == chosen_split_rank)
					fprintf(stderr,	"rank %d :::::[MPISORT][BITONIC 2] time spent = %f s\n",
											split_rank, MPI_Wtime() - time_count);
				size_t k1;
				size_t tmp2 = 0;
				for (k1 = 1; k1 < max_num_read; k1++){
					assert(local_reads_coordinates_sorted[k1-1] <= local_reads_coordinates_sorted[k1]);
					local_dest_rank_sorted[k1]= split_rank;
				}
				/*
				for (k1 = 0; k1 < max_num_read; k1++){
					fprintf(stderr,	"rank %d :::::[MPISORT][BITONIC 2]  local_reads_coordinates_sorted[%zu]= %zu s\n",
											split_rank, k1, local_reads_coordinates_sorted[k1]);

					fprintf(stderr,	"rank %d :::::[MPISORT][BITONIC 2]  local_source_rank_sorted[%zu]= %d s\n",
											split_rank, k1, local_source_rank_sorted[k1]);							
				}
				*/
				size_t *local_offset_dest_sorted = malloc(max_num_read*sizeof(size_t));
				size_t last_local_offset = 0;


				// We compute the local_dest_offsets_sorted
				size_t local_total_offset = 0;

				for (k1 = 0; k1 <  max_num_read; k1++){
					local_offset_dest_sorted[k1] = local_reads_sizes_sorted[k1];
					local_total_offset += local_reads_sizes_sorted[k1];

				}

				//we make the cumulative sum of all offsets
				for (k1 = 1; k1 < max_num_read; k1++){
					local_offset_dest_sorted[k1] = local_offset_dest_sorted[k1 - 1] + local_offset_dest_sorted[k1];
				}

				//we exchange the last destination offset
				last_local_offset = local_offset_dest_sorted[max_num_read-1];


				//number of block to send
				int blocksize = 1;

				MPI_Offset *y  = calloc(split_size, sizeof(MPI_Offset));
				MPI_Offset *y2 = calloc(split_size + 1, sizeof(MPI_Offset));

				//we wait all processors

				MPI_Gather(&last_local_offset, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, split_comm);

				if (split_rank ==0){
					for (k1 = 1; k1 < (split_size + 1); k1++) {
						y2[k1] = y[k1-1];
					}
				}

				if (split_rank ==0){
					for (k1 = 1; k1 < (split_size +1); k1++) {
						y2[k1] = y2[k1-1] + y2[k1];
					}
				}

				size_t offset_to_add = 0;
				MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &offset_to_add, 1, MPI_LONG_LONG_INT, 0, split_comm);

				free(y);
				free(y2);

				//we add offset of the previous rank
				for (k1 = 0; k1 < max_num_read; k1++){
					if (local_reads_sizes_sorted[k1] != 0)
						local_offset_dest_sorted[k1] += offset_to_add;
						//local_offset_dest_sorted[k1]=max_num_read*split_rank+k1+1;
					else
						local_offset_dest_sorted[k1] = 0;
				}


				/*
				for (k1 = 0; k1 < max_num_read; k1++){

					fprintf(stderr, "\n");

					fprintf(stderr,	"rank %d :::::[MPISORT][BITONIC 2]  local_reads_coordinates_sorted[%zu]= %zu s\n",
											split_rank, k1, local_reads_coordinates_sorted[k1]);

					fprintf(stderr,	"rank %d :::::[MPISORT][BITONIC 2]  local_source_rank_sorted[%zu]= %d s\n",
											split_rank, k1, local_source_rank_sorted[k1]);							
				

					fprintf(stderr,	"rank %d :::::[MPISORT][BITONIC 2]  local_offset_dest_sorted[%zu]= %d s\n",
											split_rank, k1, local_offset_dest_sorted[k1]);							
				
					fprintf(stderr, "\n");
				}
				*/

				/*
				 * we update destination rank according to
				 * original number of reads read.
				 *
				 */

				//we compute the new rank dest according to max_num_read
				size_t previous_num_reads_per_job[dimensions];
				//we create a vector of size split_size with previous reads per job
				MPI_Allgather(&first_local_readNum , 1, MPI_LONG_LONG_INT, previous_num_reads_per_job , 1, MPI_LONG_LONG_INT, split_comm);

				// we compute the position of of the read in the first
				// reference without the zero padding of bitonic
				size_t pos_ref0 = 0;

				//we need the number of zeros we add for the padding
				size_t N0 = max_num_read*dimensions - total_num_read;

				int new_rank = 0;
				int previous_rank = 0;
				// we compute the new rank for
				// the reads sorted by offset destination
				size_t h = 0;


				pos_ref0 = max_num_read*split_rank - N0;
				for(j = 0; j < max_num_read; j++) {
					if ( local_reads_sizes_sorted[j] != 0){
						int new_rank = chosen_split_rank;
						pos_ref0 = (max_num_read*split_rank +j) - N0;
						if (pos_ref0 >= 0) {
							size_t tmp2 = 0;
							for (h = 0; h < dimensions; h++){
								tmp2 += previous_num_reads_per_job[h];
								if ( pos_ref0 < tmp2)  {
									new_rank = h;
									break;
									}
								}
							previous_rank = local_dest_rank_sorted[j];
							local_dest_rank_sorted[j] = new_rank;
						}
					}
				}

				MPI_Barrier(split_comm);

				size_t offset  = 0;
				size_t numItems = 0;
				size_t num_read_for_bruck = 0;
				int *p = local_reads_sizes_sorted;
				if (p[0] != 0) {offset = 0;};
				if (p[max_num_read -1] == 0){offset = max_num_read;}
				else {while ((*p == 0) && (offset < max_num_read )){ offset++; p++;}}

				/*
				 * REMOVE ZERO PADDING BEFORE BRUCK
				 *
				 */

				time_count = MPI_Wtime();

				if (offset > 0){

					// we remove zeros in the vector we have 2 cases
					// the first offset <  max_num_read
					// and the entire vector is null
					if ( offset < max_num_read ){

						numItems = max_num_read - offset;

						local_reads_coordinates_sorted_trimmed_for_bruck    = malloc(numItems * sizeof(size_t));
						local_offset_source_sorted_trimmed_for_bruck        = malloc(numItems * sizeof(size_t));
						local_offset_dest_sorted_trimmed_for_bruck			= malloc(numItems * sizeof(size_t));
						local_reads_sizes_sorted_trimmed_for_bruck          = malloc(numItems * sizeof(int));
						local_dest_rank_sorted_trimmed_for_bruck            = malloc(numItems * sizeof(int));
						local_source_rank_sorted_trimmed_for_bruck 		    = malloc(numItems * sizeof(int));
						size_t y=0;

						for (y = 0; y < numItems; y++){

							local_reads_coordinates_sorted_trimmed_for_bruck[y]    = local_reads_coordinates_sorted[y+offset];
							local_offset_source_sorted_trimmed_for_bruck[y]        = local_offset_source_sorted[y+offset];
							local_offset_dest_sorted_trimmed_for_bruck[y]		   = local_offset_dest_sorted[y+offset];
							local_reads_sizes_sorted_trimmed_for_bruck[y]          = local_reads_sizes_sorted[y+offset];
							local_dest_rank_sorted_trimmed_for_bruck[y]            = local_dest_rank_sorted[y+offset];
							local_source_rank_sorted_trimmed_for_bruck[y] 		   = local_source_rank_sorted[y+offset];
						}

						num_read_for_bruck = numItems;

						/*
						 *
						 * FOR DEBUG
						 

						for(y = 0; y < num_read_for_bruck; y++){
							assert( local_reads_sizes_sorted_trimmed_for_bruck[y] 		!= 0 );
							assert( local_source_rank_sorted_trimmed_for_bruck[y] 		< dimensions);
							assert( local_dest_rank_sorted_trimmed_for_bruck[y]   		< dimensions);
							assert( local_offset_source_sorted_trimmed_for_bruck[y] 	!= 0);
							assert( local_offset_dest_sorted_trimmed_for_bruck[y] 	    != 0);
							assert( local_reads_coordinates_sorted_trimmed_for_bruck[y] != 0);
						}
						*/

					}
					else{

						numItems = 0;
						local_reads_coordinates_sorted_trimmed_for_bruck    = malloc(numItems * sizeof(size_t));
						local_offset_source_sorted_trimmed_for_bruck        = malloc(numItems * sizeof(size_t));
						local_offset_dest_sorted_trimmed_for_bruck          = malloc(numItems * sizeof(size_t));
						local_reads_sizes_sorted_trimmed_for_bruck          = malloc(numItems * sizeof(int));
						local_dest_rank_sorted_trimmed_for_bruck            = malloc(numItems * sizeof(int));
						local_source_rank_sorted_trimmed_for_bruck 		    = malloc(numItems * sizeof(int));
						num_read_for_bruck = 0;
					}
				}
				else {

					numItems = local_readNum;
					local_reads_coordinates_sorted_trimmed_for_bruck    = malloc(local_readNum * sizeof(size_t));
					local_offset_source_sorted_trimmed_for_bruck        = malloc(local_readNum * sizeof(size_t));
					local_offset_dest_sorted_trimmed_for_bruck          = malloc(local_readNum * sizeof(size_t));
					local_reads_sizes_sorted_trimmed_for_bruck          = malloc(local_readNum * sizeof(int));
					local_dest_rank_sorted_trimmed_for_bruck            = malloc(local_readNum * sizeof(int));
					local_source_rank_sorted_trimmed_for_bruck 		    = malloc(local_readNum * sizeof(int));

					size_t y=0;
					for (y = 0; y < local_readNum; y++){

						local_reads_coordinates_sorted_trimmed_for_bruck[y]    = local_reads_coordinates_sorted[y];
						local_offset_source_sorted_trimmed_for_bruck[y]        = local_offset_source_sorted[y];
						local_offset_dest_sorted_trimmed_for_bruck[y]          = local_offset_dest_sorted[y];
						local_reads_sizes_sorted_trimmed_for_bruck[y]          = local_reads_sizes_sorted[y];
						local_dest_rank_sorted_trimmed_for_bruck[y]            = local_dest_rank_sorted[y];
						local_source_rank_sorted_trimmed_for_bruck[y] 		   = local_source_rank_sorted[y];
					}

					num_read_for_bruck = numItems;

					/*
					 *
					 * FOR DEBUG
					 
					for(y = 0; y < num_read_for_bruck; y++){
						assert( local_reads_sizes_sorted_trimmed_for_bruck[y] 		!= 0 );
						assert( local_source_rank_sorted_trimmed_for_bruck[y] 		< dimensions);
						assert( local_dest_rank_sorted_trimmed_for_bruck[y]   		< dimensions);
						assert( local_offset_source_sorted_trimmed_for_bruck[y] 	!= 0);
						assert( local_offset_dest_sorted_trimmed_for_bruck[y] 	    != 0);
						assert( local_reads_coordinates_sorted_trimmed_for_bruck[y] != 0);
					}
					*/
					
				}

				free(local_reads_coordinates_sorted);
				free(local_offset_source_sorted);
				free(local_offset_dest_sorted);
				free(local_reads_sizes_sorted);
				free(local_dest_rank_sorted);
				free(local_source_rank_sorted);


				if (split_rank == chosen_split_rank)
					fprintf(stderr,	"rank %d :::::[MPISORT][TRIMMING] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);

				/*
				 * We do a Bruck on rank of origin reading
				 */

				size_t m=0;
				int num_proc = dimensions;
				size_t *number_of_reads_by_procs = calloc( dimensions, sizeof(size_t));

				//fprintf(stderr,	"rank %d :::::[MPISORT] num_read_for_bruck = %zu \n", split_rank, num_read_for_bruck);

				for(m = 0; m < num_read_for_bruck; m++){
					 //assert(new_pbs_orig_rank_off_phase1[m] < dimensions);
					 //assert(new_pbs_dest_rank_phase1[m] < dimensions);
					 number_of_reads_by_procs[local_source_rank_sorted_trimmed_for_bruck[m]]++;
				}

				int *local_source_rank_sorted_trimmed_for_bruckv2 = malloc( num_read_for_bruck * sizeof(int));

				for(m = 0; m < num_read_for_bruck; m++){
					local_source_rank_sorted_trimmed_for_bruckv2[m] = local_source_rank_sorted_trimmed_for_bruck[m];
				}

				size_t count6 = 0;
				for(m = 0; m < dimensions; m++){
					count6 += number_of_reads_by_procs[m];
				}

				assert( count6 == num_read_for_bruck );
				MPI_Barrier(split_comm);

				size_t **reads_coordinates 		= malloc(sizeof(size_t *) * dimensions);
				size_t **local_source_offsets 	= malloc(sizeof(size_t *) * dimensions);
				size_t **dest_offsets 			= malloc(sizeof(size_t *) * dimensions);
				int **read_size 				= malloc(sizeof(int *) * dimensions);
				int **dest_rank 				= malloc(sizeof(int *) * dimensions);
				int **source_rank				= malloc(sizeof(int *) * dimensions);

				/*
				 * We send in order
				 *
				 * local_offset_source_sorted_trimmed_for_bruck
				 * local_dest_rank_sorted_trimmed_for_bruck
				 * local_reads_coordinates_sorted_trimmed_for_bruck
				 * local_reads_sizes_sorted_trimmed_for_bruck
				 *
				 */

				COMM_WORLD = split_comm;
				time_count = MPI_Wtime();

				bruckWrite3(split_rank,
							dimensions,
							count6,
							number_of_reads_by_procs,
							local_source_rank_sorted_trimmed_for_bruckv2,
							local_offset_source_sorted_trimmed_for_bruck,     //offset sources
							&local_source_offsets,
							local_dest_rank_sorted_trimmed_for_bruck,     	  //destination rank
							&dest_rank,
							local_reads_coordinates_sorted_trimmed_for_bruck, //reads coordinates
							&reads_coordinates,
							local_reads_sizes_sorted_trimmed_for_bruck,       //read size
							&read_size,
							local_source_rank_sorted_trimmed_for_bruck,		  //source rank
							&source_rank,
							local_offset_dest_sorted_trimmed_for_bruck,
							&dest_offsets
				);

				if (split_rank == chosen_split_rank)
					fprintf(stderr,	"rank %d :::::[MPISORT][BRUCK 3] time spent = %f s\n",
							split_rank, MPI_Wtime() - time_count);


				time_count = MPI_Wtime();

				free(local_reads_coordinates_sorted_trimmed_for_bruck);
				free(local_dest_rank_sorted_trimmed_for_bruck);
				free(local_reads_sizes_sorted_trimmed_for_bruck);
				free(local_offset_source_sorted_trimmed_for_bruck);
				free(local_offset_dest_sorted_trimmed_for_bruck);
				free(local_source_rank_sorted_trimmed_for_bruck);
				free(local_source_rank_sorted_trimmed_for_bruckv2);

				local_reads_coordinates_sorted_trimmed 	  = malloc(first_local_readNum * sizeof(size_t));
				local_offset_source_sorted_trimmed   	  = malloc(first_local_readNum * sizeof(size_t));
				local_offset_dest_sorted_trimmed   	  	  = malloc(first_local_readNum * sizeof(size_t));
				local_dest_rank_sorted_trimmed   		  = malloc(first_local_readNum * sizeof(int));
				local_source_rank_sorted_trimmed		  = malloc(first_local_readNum * sizeof(int));
				local_reads_sizes_sorted_trimmed		  = malloc(first_local_readNum * sizeof(int));

				if (split_rank == chosen_split_rank)
					fprintf(stderr,	"rank %d :::::[MPISORT][FREE + MALLOC] time spent = %f s\n",
											split_rank, MPI_Wtime() - time_count);
				/*
				 * GET DATA AFTER BRUCK
				 *
				 */

				j=0;
				size_t k = 0;

				for(m = 0; m < num_proc; m++)
				{
					for(k = 0; k < number_of_reads_by_procs[m]; k++)
					{
						
						local_offset_dest_sorted_trimmed[k + j] 		= dest_offsets[m][k];
						local_dest_rank_sorted_trimmed[k + j] 			= dest_rank[m][k];
						local_reads_sizes_sorted_trimmed[k + j] 		= read_size[m][k];
						local_offset_source_sorted_trimmed[k + j] 		= local_source_offsets[m][k];
						local_reads_coordinates_sorted_trimmed[k + j] 	= reads_coordinates[m][k];
						local_source_rank_sorted_trimmed[k + j] 		= source_rank[m][k];

					}
					free(dest_offsets[m]);
					free(dest_rank[m]);
					free(read_size[m]);
					free(local_source_offsets[m]);
					free(reads_coordinates[m]);
					free(source_rank[m]);
					j += number_of_reads_by_procs[m];
				}


				free(number_of_reads_by_procs);
				if (dest_rank != NULL)
					free(dest_rank);
				if (read_size != NULL)
					free(read_size);
				if (local_source_offsets != NULL)
					free(local_source_offsets);
				if (reads_coordinates != NULL)
					free(reads_coordinates);
				if (source_rank != NULL)
					free(source_rank);
				if (dest_offsets != NULL)
					free(dest_offsets);

				local_readNum = first_local_readNum;


				/*
				 *
				 * FOR DEBUG
				 
				for ( j = 0; j < local_readNum; j++){
					assert ( local_reads_coordinates_sorted_trimmed[j]    != 0 );
					assert ( local_offset_source_sorted_trimmed[j]        != 0 );
					assert ( local_offset_dest_sorted_trimmed[j]   		  != 0 );
					assert ( local_reads_sizes_sorted_trimmed 			  != 0 );
					assert ( local_dest_rank_sorted_trimmed[j]            < split_size );
					assert ( local_source_rank_sorted_trimmed[j] 		  < split_size );
				}
				*/

				free(local_reads_coordinates_sorted_trimmed);

				if (split_rank == chosen_split_rank)
					fprintf(stderr,	"rank %d :::::[MPISORT] we call write SAM \n", split_rank);

				malloc_trim(0);

				time_count = MPI_Wtime();

				writeSam(
					split_rank,
					output_dir,
					header,
					local_readNum,
					total_reads_by_chr,
					chrNames[i],
					reads[i],
					split_size,
					split_comm,
					chosen_split_rank,
					file_name,
					mpi_file_split_comm,
					finfo,
					compression_level,
					local_offset_dest_sorted_trimmed,
					local_offset_source_sorted_trimmed,
					local_reads_sizes_sorted_trimmed,
					local_dest_rank_sorted_trimmed,
					local_source_rank_sorted_trimmed,
					local_data,
					goff[rank],
					first_local_readNum,
					uniq_chr,
					write_format,
					merge,
					file_name_merge
				);

				if (split_rank == chosen_split_rank){
					fprintf(stderr,	"rank %d :::::[MPISORT][WRITESAM] chromosom %s :::  %f seconds\n\n\n",
							split_rank, chrNames[i], MPI_Wtime() - time_count);

				}
			}
			else{

				/*
				 * We are in the case the number of cpu is
				 * not a power of 2
				 *
				 *
				 */

				parallel_sort_any_dim(
						dimensions, 				//dimension for parabitonic
						local_readNum,
						split_rank,
						split_size,
						reads,
						i, 							//chromosom number
						chosen_split_rank,
						split_comm,
						localReadNumberByChr,
						local_data,
						file_name,
						output_dir,
						finfo,
						compression_level,
						total_reads_by_chr,
						goff[rank],
						headerSize,
						header,
						chrNames[i],
						mpi_file_split_comm,
						uniq_chr,
						write_format,
                				merge,
                				file_name_merge
					);

			} //end if dimensions < split_rank

		} //if ((local_color == 0) && (i < (nbchr - 2))) //in the splitted dimension
		else{
			//we do nothing here
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

	}// end loop upon chromosoms (line 665)


	free(goff);
	if (!uniq_chr)	free(local_data);

	MPI_Barrier(MPI_COMM_WORLD);

	free(header); //ok
	free(localReadNumberByChr); //ok

	for(i = 0; i < nbchr; i++){
		free(chrNames[i]);
	}

	free(reads); //ok
	free(readNumberByChr); //ok

	malloc_trim(0);
	res = MPI_Finalize();
	assert(res == MPI_SUCCESS);

	return 0;
}



static void usage(const char *prg) {

	fprintf(stderr, "program: %s is a MPI version for sorting SAM file\n"
		"version: %s\n"
		"\nusage : mpirun -n TOTAL_PROC %s SAM_FILE OUTPUT_DIRECTORY -q QUALITY -n \n"
        "\n\tTOTAL_PROC tells how many cores will be used by MPI to parallelize the computation.\n"
        "\noptions:\n"
        "\n\t-p if the read are paired-end (by defaut reads are single-end)\n"
	"\n\t-u if the file contains only one chromosome for instance results from mpiBwaByChr (by defaut all chromosomes are present)\n"
        "\n\t-q INTEGER\n"
        "\t     filters the reads according to their quality. Reads quality under the\n"
        "\t     threshold are ignored in the sorting results. Default is 0 (all reads are kept).\n"
        "\n\t-n\n"
        "\t     sorts the read by their query name.\n"
	"\n\t-s\n"
	"\t	write the output in sam format.\n"
        "\ninput: input file is a sam file of paired or single reads\n"
        "\noutput: set of gz files with\n"
        "\t* one per chromosome (e.g. chr11.gz)\n"
        "\t* one for discordant reads (discordant.gz): discordants reads are reads \n"
        "\t  where one pair aligns on a chromosome and the other pair aligns on \n"
        "\t  another chromosome \n"
        "\t* one for unmapped reads (unmapped.gz): unmapped reads are reads without \n"
        "\t  coordinates on any chromosome \n"
		"\nexample : mpirun -n 4 %s  HCC1187C_70K_READS.sam ${HOME}/mpiSORTExample -q 0 -n \n"
        "\nFor more detailed documentation visit:\n"
        "\thttps://github.com/bioinfo-pf-curie/mpiSORT\n"
        "\nCopyright (C) 2020  Institut Curie <http://www.curie.fr> \n"
        "\nThis program comes with ABSOLUTELY NO WARRANTY. \n"
        "This is free software, and you are welcome to redistribute it \n"
        "under the terms of the CeCILL License. \n"
		"\ncontact: Frederic Jarlier (frederic.jarlier@curie.fr) \n"
		,prg, VERSION, prg, prg);

	return; }
