/*
   This file is part of mpiSORT
   
   Copyright Institut Curie 2022
   
   This software is a computer program whose purpose is to sort SAM file.
   
   You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
   
   The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
   
   The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/

/*
   Module:
     parser.c

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

#include "parser.h"

char parse_mode;

size_t hash_name(char *begin_name, char *end_name, int max)
 {
 	int i, j;
 	size_t result = 0;
	size_t name_length = end_name - begin_name;
	char *line = begin_name; 

 	for(i = 0, j = 0; i < name_length; i++) {
 		if( (line[i] >= 'A' && line[i] <= 'Z') || (line[i] >= 'a' && line[i] <= 'z') )
 			j = i;
 	}

 	i = j + 1;
 	for(j = 0; i < name_length && j < max; i++) {
 		if(line[i] >= '0' && line[i] <= '9') {
 			char tmp;
 			tmp = line[i] - '0';
 			result <<= 4;
 			result += tmp;
 			j++;
 		}
 	}
 	return result;
 }

/*
size_t hash_name(char *s, int max){
	khint_t h = (khint_t)*s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + (khint_t)*s;
	return h;
}
*/

void init_goff(MPI_File mpi_filed,unsigned int headerSize,size_t fsize,int numproc,int rank, size_t *goff){


	char * current_line = NULL;
	MPI_Status status;
	int i = 0;
	int j = 0;

	size_t lsize = fsize/numproc;
	goff[0]=headerSize;
	for(i=1;i<numproc;i++){
		goff[i]=lsize*i+headerSize;
	}

	goff[numproc]=fsize;

	for(i=1;i<numproc;i++)
	{
		current_line =(char*)calloc(1024*1024,sizeof(char));
		MPI_File_read_at(mpi_filed, (MPI_Offset)goff[i], current_line, 1024*1024, MPI_CHAR, &status);
		assert(strlen(current_line) != 0);
		j=0;
		while(j<fsize && current_line[j] != '\n'){j++;}
		goff[i]+=(j+1);
		free(current_line);
	}


}



int parser_single(char *localData, size_t size2compute, int rank, size_t start_offset, unsigned char threshold,
		int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads , CHTbl *ref_htbl){

		#define _read_token_tab(_p) (_p); do { char *tab = strchr((_p), '\t'); if (!tab) goto err_ret;  (_p) = tab + 1; } while (0)
                #define _read_token_ret(_p) (_p); do { char *tab = strchr((_p), '\n'); if (!tab) goto err_ret;  (_p) = tab + 1; } while (0)

                char *currentCarac;
		char *p = localData;
                char *p1 = localData;
                char *q;
                char *q1, *q2;
                char *name;
                char *flag;
                char *chr_txt;
                char *mchr_txt;
                char *coord_txt;
                char *qual_txt;
		char chr_name[100];
		//for debug
                //char currentLine[MAX_LINE_SIZE];
                unsigned char quality;
                unsigned int i, chr, nbchr = 0;
                int lastChr = -1;
                int next;
		size_t total_computed = 0;
                size_t lineSize, offset_read_in_source_file;
                size_t coord;
                size_t *readNumberByChr;
                size_t counter = 0;
                Read **reads = *preads;

		//for debug
                //for(i=0;i<MAX_LINE_SIZE;i++) currentLine[i]=0;
                
		 //we take the first line *
                 //before calling parsepaired, we know that localdata is at the begining of a read
                 //next = tokenizer(localData,'\n', currentLine);
                 offset_read_in_source_file = start_offset;
                 nbchr = nbchrom;
                 readNumberByChr = (size_t*)calloc(nbchr, sizeof(size_t));

		 while( total_computed < size2compute) {
		
			for (i = 0; i<100; i++) chr_name[i]=0;
                        q = _read_token_ret(p);
			
			//for debug
	                //memcpy(currentLine, q , p -q );
                        //currentLine[p-q]=0;
                        assert(*p != '\n');
                        total_computed += p - q;
                        lineSize = p - q;
 			//get the name
                        name = _read_token_tab(q);
			if (parse_mode == MODE_NAME) { q1 = name; q2 = q; coord = hash_name(q1, q2, 16);}
                        //get the flag
                        flag = _read_token_tab(q);

                        //get the chr name)
                        chr_txt = _read_token_tab(q);
                        if ( *chr_txt == '*' ) chr = (nbchr-1);
                        else{
                                memcpy(chr_name, chr_txt , q - chr_txt - 1);
                                //look the index of the chromosom in hash table
                                chr = chtbl_lookup(ref_htbl, chr_name);
                                assert(chr >= 0);
                        }

			//GO TO COORD
                        if (parse_mode != MODE_NAME) {
                                coord_txt = _read_token_tab(q);
                                coord = strtoull(coord_txt, &q, 10);
                        }

                        qual_txt = _read_token_tab(q);
                        quality = strtoull(q, &qual_txt, 10);

			//fprintf(stderr, "currentLine = %s \n chr_name = %s :::: chr = %d :::: coord = %zu  \n ", currentLine, chr_name, chr, coord);

			if (chr < nbchr-2){                                
				 if(quality >= threshold){

                                       reads[chr]->next = malloc(sizeof(Read));
                                       reads[chr]->next->coord = coord;
                                       reads[chr]->next->quality = quality;
                                       reads[chr]->next->offset_source_file=offset_read_in_source_file;
                                       reads[chr]->next->offset = lineSize;
                                       reads[chr] = reads[chr]->next;
                                       readNumberByChr[chr]++;
                                }

			}
			else {
                                reads[nbchr-1]->next = malloc(sizeof(Read));
                                reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
                                reads[nbchr-1]->next->offset = lineSize;
                                reads[nbchr-1] = reads[nbchr-1]->next;
                                readNumberByChr[nbchr-1]++;
                        }
                                 			
			offset_read_in_source_file += lineSize;
                        q = p;
                        counter++;
                       
                }

                fprintf(stderr, "rank %d ::: counter = %zu \n", rank, counter);

                for(i=0;i<nbchr;i++){
                        preadNumberByChr[0][i] += readNumberByChr[i];
                }
                free(readNumberByChr);

		return 0;

                err_ret:
                        return -2;

}


int parser_paired_uniq(char *localData, size_t size2compute, int rank, size_t start_offset, unsigned char threshold,
                int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads, int print){


		#define _read_token_tab(_p) (_p); do { char *tab = strchr((_p), '\t'); if (!tab) goto err_ret;  (_p) = tab + 1; } while (0)
		#define _read_token_ret(_p) (_p); do { char *tab = strchr((_p), '\n'); if (!tab) goto err_ret;  (_p) = tab + 1; } while (0) 

		char *currentCarac;
		char *p = localData;
		char *p1 = localData;
		char *q;
		char *q1, *q2;
		char *name;
		char *flag;
		char *chr_txt;
		char *mchr_txt; 
		char *coord_txt;
		char *qual_txt;
                unsigned char quality;
                unsigned int i, chr, nbchr = 0, mchr;
                int lastChr = -1;
                int next;
		size_t total_computed = 0;
                size_t lineSize, offset_read_in_source_file;
                size_t coord;
                size_t *readNumberByChr;
                size_t counter = 0;
                Read **reads = *preads;

		 
		 offset_read_in_source_file = start_offset;
		 nbchr = nbchrom;
		 readNumberByChr = (size_t*)calloc(nbchr, sizeof(size_t));
		 
		 while( total_computed < size2compute){

			q = _read_token_ret(p);
			assert(*p != '\n');	
                 	total_computed += p - q;
                 	lineSize = p - q;

			//we pass the name 
			name = _read_token_tab(q);
			if (parse_mode == MODE_NAME) { q1 = name; q2 = q; coord = hash_name(q1, q2, 16);}
			//pass the flag
			flag = _read_token_tab(q);
			//pass the chr
			chr_txt = _read_token_tab(q);  
			chr = 0;                        
			
			if (parse_mode != MODE_NAME) {
                		coord_txt = _read_token_tab(q);
                                coord = strtoull(coord_txt, &q, 10);
                                
                        }
                        
			qual_txt = _read_token_tab(q);
                        quality = strtoull(q, &qual_txt, 10);
			if (parse_mode != MODE_NAME) _read_token_tab(q);		
			
			mchr_txt = _read_token_tab(q);
                        if( *q == '=') mchr = chr;
                        else if( *q == '*') mchr = (nbchr-1);
                        else mchr = 0;

        		if (chr < nbchr-2){
				if(quality >= threshold){

                                                reads[chr]->next = malloc(sizeof(Read));
                                                reads[chr]->next->coord = coord;
                                                reads[chr]->next->quality = quality;
                                                reads[chr]->next->offset_source_file=offset_read_in_source_file;
                                                reads[chr]->next->offset = lineSize;
                                                reads[chr] = reads[chr]->next;
                                                readNumberByChr[chr]++;
                                        }
			}
			else {
				//we found unmapped pairs reads
				reads[nbchr-1]->next = malloc(sizeof(Read));
				reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
				reads[nbchr-1]->next->offset = lineSize;
				reads[nbchr-1] = reads[nbchr-1]->next;
				readNumberByChr[nbchr-1]++;
				
			}
			
			offset_read_in_source_file += lineSize;
                        q = p;			
			counter++;
		}

		fprintf(stderr," rank = %d ::::counter = %zu \n", rank, counter);

		//fprintf(stderr, "rank %d ::: read number = %zu \n", rank, counter);	
                //fprintf(stderr, "rank %d ::: total_computed = %zu \n", rank, total_computed);
		//fprintf(stderr, "rank %d ::: size2compute = %zu \n", rank, size2compute);
		assert(total_computed == size2compute);
                for(i=0;i<nbchr;i++){
                        preadNumberByChr[0][i] += readNumberByChr[i];
                }
                free(readNumberByChr);
		return 0;
		err_ret:
			return -2;
}



int parser_paired(char *localData, size_t size2compute, int rank, size_t start_offset, unsigned char threshold,
		int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads, CHTbl *ref_htbl){

		#define _read_token_tab(_p) (_p); do { char *tab = strchr((_p), '\t'); if (!tab) goto err_ret;  (_p) = tab + 1; } while (0)
                #define _read_token_ret(_p) (_p); do { char *tab = strchr((_p), '\n'); if (!tab) goto err_ret;  (_p) = tab + 1; } while (0) 

		char *currentCarac;
		char *p = localData;
                char *p1 = localData;
                char *q;
                char *q1;
		char *q2;
                char *name;
                char *flag;
                char *chr_txt;
		char *mchr_txt;
                char *coord_txt;
                char *qual_txt;
		char *cigar, cigar_str[100];
		char chr_name[100];
		char mchr_name[100];
		unsigned char quality;
		unsigned int i, chr, nbchr = 0, mchr;
		int lastChr = -1;
		int next;
		size_t lineSize, offset_read_in_source_file;
		size_t coord;
		size_t *readNumberByChr;
		size_t counter = 0;
		size_t total_computed = 0;
		Read **reads = *preads;

		offset_read_in_source_file = start_offset;
		nbchr = nbchrom;
		readNumberByChr = (size_t*)calloc(nbchr, sizeof(size_t));

		while( total_computed < size2compute ){

			for (i = 0; i<100; i++){chr_name[i]=0;mchr_name[i]=0;}
			q = _read_token_ret(p);

			assert(*p != '\n');
                        total_computed += p - q;
			lineSize = p - q;
			//get the name
			name = _read_token_tab(q);
			if (parse_mode == MODE_NAME) { q1 = name; q2 = q; coord = hash_name(q1, q2, 16);}
                        
			//get the flag
			flag = _read_token_tab(q);
			//get the chr name)
			chr_txt = _read_token_tab(q);
			if ( *chr_txt == '*' ) chr = (nbchr-1);
			else{
				memcpy(chr_name, chr_txt , q - chr_txt - 1);
				//look the index of the chromosom in hash table
				chr = chtbl_lookup(ref_htbl, chr_name);
				assert(chr >= 0);   				 
			}
			
			//GO TO COORD
			//coord_txt = _read_token_tab(q);
			if (parse_mode != MODE_NAME) {
				coord_txt = _read_token_tab(q);
				coord = strtoull(coord_txt, &q, 10);
			}
			//got to qual 
			qual_txt = _read_token_tab(q);
			quality = strtoull(q, &qual_txt, 10);
			//goto cigar
			if (parse_mode != MODE_NAME) _read_token_tab(q);
			cigar =  _read_token_tab(q);
			//goto mate mchr
			mchr_txt = _read_token_tab(q);

                        if( *mchr_txt == '=') mchr = chr;
                        else if( *mchr_txt == '*') mchr = (nbchr-1);
			else {
				memcpy(mchr_name, mchr_txt , q - mchr_txt - 1);
				mchr = chtbl_lookup(ref_htbl, mchr_name);
			}
			//first we check if reads mapped on the same chromosome
			if ((chr < nbchr-2) && (chr == mchr)){
					//then we found concordant reads
					if(quality >= threshold){

						reads[chr]->next = malloc(sizeof(Read));
						reads[chr]->next->coord = coord;
						reads[chr]->next->quality = quality;
						reads[chr]->next->offset_source_file=offset_read_in_source_file;
						reads[chr]->next->offset = lineSize;
						reads[chr] = reads[chr]->next;
						readNumberByChr[chr]++;
					}
			}
			else if ((chr < (nbchr-2)) && ( mchr < (nbchr -2)) && ( chr != mchr)){

					//we found discordant reads
					if(quality >= threshold){

                                                reads[chr]->next = malloc(sizeof(Read));
                                                reads[chr]->next->coord = coord;
                                                reads[chr]->next->quality = quality;
                                                reads[chr]->next->offset_source_file=offset_read_in_source_file;
                                                reads[chr]->next->offset = lineSize;
                                                reads[chr] = reads[chr]->next;
                                                readNumberByChr[chr]++;
                                        }

			}
			else if ((chr == (nbchr-1)) && ( mchr < (nbchr -2))){

					//we found discordant reads with one pair unmapped
					if(quality >= threshold){

                                                reads[mchr]->next = malloc(sizeof(Read));
                                                reads[mchr]->next->coord = coord;
                                                reads[mchr]->next->quality = quality;
                                                reads[mchr]->next->offset_source_file=offset_read_in_source_file;
                                                reads[mchr]->next->offset = lineSize;
                                                reads[mchr] = reads[mchr]->next;
                                                readNumberByChr[mchr]++;
                                        }
			}
			else if ((mchr == (nbchr-1)) && ( chr < (nbchr -2))){

					//we found discordant reads with one pair unmapped
					
					if(quality >= threshold){

                                                reads[chr]->next = malloc(sizeof(Read));
                                                reads[chr]->next->coord = coord;
                                                reads[chr]->next->quality = quality;
                                                reads[chr]->next->offset_source_file=offset_read_in_source_file;
                                                reads[chr]->next->offset = lineSize;
                                                reads[chr] = reads[chr]->next;
                                                readNumberByChr[chr]++;
                                        }


			}
			else{
					//we found unmapped pairs reads
					//fprintf(stderr, "\n \n");
					//fprintf(stderr, "rank %d unmapped currentLine =  %s \n", rank, currentLine);
					//fprintf(stderr, "rank %d chr =  %d \n", rank, chr);
					//fprintf(stderr, "rank %d mchr =  %d \n", rank, mchr);
					//fprintf(stderr, "rank %d quality =  %s \n", rank, quality);
					reads[nbchr-1]->next = malloc(sizeof(Read));
					reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
					reads[nbchr-1]->next->offset = lineSize;
					reads[nbchr-1] = reads[nbchr-1]->next;
					readNumberByChr[nbchr-1]++;
			}


			offset_read_in_source_file += lineSize;
                        q = p;
                        counter++;
			//we update the offset_read_in_source_file
			/*for(i=0;i<MAX_LINE_SIZE;i++){
				currentLine[i]=0;
			}
			next = tokenizer(NULL, '\n', currentLine);
			*/	
		}

		fprintf(stderr, "rank %d ::: counter = %zu \n", rank, counter);
		assert(total_computed == size2compute);
		for(i=0;i<nbchr;i++){
			preadNumberByChr[0][i] += readNumberByChr[i];
		}
		free(readNumberByChr);
		return 0;

		err_ret:
                	return -2;


}

int clear_htable(CHTbl *ref_htbl){
	chtbl_destroy(ref_htbl);
}


int getChr(char *str, char** chrNames, int nbchr){
	int i=0, found=0, size;
	char *str1 = str, *str2;

	str2 = str1+1;

	for(; *str2 != '\t'; str2++);

	size=strlen(str1)-strlen(str2);

	char *tmp_chr =(char*)malloc(sizeof(char)*size+1);
	tmp_chr[0]=0;

	for(i=0;i<size;i++){
		tmp_chr[i]=str1[i+1];
	}
	tmp_chr[size-1]=0;

	assert(strlen(tmp_chr) != 0);
	for(i = 0, found = 0; i < nbchr && !found; i++){
		found = !strcmp(tmp_chr, chrNames[i]);
	}

	free(tmp_chr);
	if (found) return i - 1;
	else return i = 65535;
}

void get_coordinates_and_offset_source_and_size_and_free_reads(int rank, int *local_read_rank, size_t *coordinates,
                size_t* offset, int* size, Read* data_chr, int local_readNum){

        size_t j;
        Read* chr = data_chr;
        Read* to_free = chr;
        for(j = 0; j < local_readNum; j++){
                coordinates[j] = 0;
                size[j] = 0;
                offset[j] = 0;
                local_read_rank[j]=rank;
        }
        for(j = 0; j < local_readNum; j++){
                coordinates[j] = chr->coord;
                offset[j] = chr->offset_source_file;
                size[j] = (int)chr->offset;
                to_free = chr;
                chr = chr->next;
                if (to_free) free(to_free);
        }
}

size_t init_coordinates_and_size(int rank, int *local_reads_rank, size_t *local_reads_index,
                size_t* coordinates, int* size, Read* data_chr, int local_readNum)
{
        size_t dataSize = 0;
        size_t j;
        Read* chr = data_chr;
        for(j = 0; j < local_readNum; j++){
                size[j] = 0;
                coordinates[j] = 0;
        }
        for(j = 0; j < local_readNum; j++){
                coordinates[j] = chr->coord;
                size[j] = (int)chr->offset; 
                local_reads_rank[j] = rank;
                local_reads_index[j] = j;
                dataSize += chr->offset;
                chr = chr->next;
        }
        return dataSize;
}
size_t init_coordinates_and_size2(int rank, int *local_reads_rank,
                size_t* coordinates, int* size, Read* data_chr, int local_readNum)
{
        size_t dataSize = 0;
        size_t j;
        Read* chr = data_chr;
        for(j = 0; j < local_readNum; j++){
                size[j] = 0;
                coordinates[j] = 0;
        }
        for(j = 0; j < local_readNum; j++){
                coordinates[j] = chr->coord;
                size[j] = (int)chr->offset;
                local_reads_rank[j] = rank;
                dataSize += chr->offset;

                chr = chr->next;
        }

        return dataSize;
}


void chosen_split_rank_gather_size_t(MPI_Comm split_comm, int rank, int num_proc, int master, size_t size, size_t *size_per_jobs,
                size_t *start_size_per_job, size_t *all_data, size_t *data, size_t start_index)
{
        MPI_Status status;
        size_t j, k;
        if (rank == master){
                size_t st = start_size_per_job[master];
                for (k = 0; k < size_per_jobs[master]; k++){
                        all_data[st] = data[k+start_index];
                        st++;
                }
                for(j = 0; j < num_proc; j++){
                        if (j != master && size_per_jobs[j] != 0){
                                size_t temp_buf[size_per_jobs[j]];
                                temp_buf[size_per_jobs[j]] = 0;
                                assert(temp_buf !=0 );
                                MPI_Recv(temp_buf, size_per_jobs[j], MPI_LONG_LONG_INT, j, 0, split_comm, &status);
                                size_t st = start_size_per_job[j];
                                for (k = 0; k < size_per_jobs[j]; k++){
                                        all_data[st] = temp_buf[k];
                                        st++;
                                }
                        }
                }
        }
        else{
                MPI_Send(data, size, MPI_LONG_LONG_INT, master,  0, split_comm);
        }
}




