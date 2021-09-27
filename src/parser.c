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

#include <unistd.h>

#include "parser.h"

char parse_mode;

size_t hash_name(char * line, int max)
 {
 	int i, j;
 	size_t result = 0;
 	for(i = 0, j = 0; i < strlen(line); i++) {
 		if( (line[i] >= 'A' && line[i] <= 'Z') || (line[i] >= 'a' && line[i] <= 'z') )
 			j = i;
 	}

 	i = j + 1;
 	for(j = 0; i < strlen(line) && j < max; i++) {
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

void init_goff(int mpi_file_descriptor, unsigned int headerSize,size_t fsize,int numproc,int rank, size_t *goff){


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
		current_line =(char*)calloc(1000,sizeof(char));

		int ret = pread(mpi_file_descriptor, current_line, 1000, goff[i]);
		assert(ret > 0);
		j=0;
		while(j<fsize && current_line[j] != '\n'){j++;}
		goff[i]+=(j+1);
		free(current_line);
	}


}



void parser_single(char *localData, int rank, size_t start_offset, unsigned char threshold,
		int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads){

                char *currentCarac;
                char currentLine[MAX_LINE_SIZE];
                unsigned char quality;
                unsigned int i, chr, nbchr = 0;
                int lastChr = -1;
                int next;
                size_t lineSize, offset_read_in_source_file;
                size_t coord;
                size_t *readNumberByChr;
                size_t counter = 0;
                Read **reads = *preads;

                for(i=0;i<MAX_LINE_SIZE;i++){
                        currentLine[i]=0;
                }
		 //we take the first line *
                 //before calling parsepaired, we know that localdata is at the begining of a read
                 next = tokenizer(localData,'\n', currentLine);
                 offset_read_in_source_file = start_offset;
                 nbchr = nbchrom;
                 readNumberByChr = (size_t*)calloc(nbchr, sizeof(size_t));

		 while(next){

                        lineSize = strlen(currentLine) + 1;
			 //we update the offset in the
			 //source file
			currentLine[lineSize - 1] = '\n';
                        currentLine[lineSize] = '\0';
			//GO TO FLAG
                        currentCarac = strstr(currentLine, "\t");
			*currentCarac = '\0';
			currentCarac++;
			 currentCarac = strstr(currentCarac+1, "\t");
                        if(currentCarac[1] == '*'){
                                chr = (nbchr-1);
                        }
                        else
                        {
                                chr = getChr(currentCarac, chrNames, nbchr);

                        }
			currentCarac = strstr(currentCarac+1, "\t");
			if (parse_mode == MODE_NAME) {
                                coord = strtoull(currentLine, NULL, strlen(currentLine));
				coord = hash_name(currentLine, 16);
				strtoull(currentCarac, &currentCarac, 10);
	                }
			else {
                                //TAKE COORD AND GO TO MAPQ
                                coord = strtoull(currentCarac, &currentCarac, 10);

                         quality = strtoull(currentCarac, &currentCarac, 10);}
			
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
			else if ((chr == 65535)){
                                reads[nbchr-2]->next = malloc(sizeof(Read));
                                reads[nbchr-2]->next->offset_source_file=offset_read_in_source_file;
                                reads[nbchr-2]->next->offset = lineSize;
                                reads[nbchr-2] = reads[nbchr-2]->next;
                                readNumberByChr[nbchr-2]++;
                        }    
			else {
                                reads[nbchr-1]->next = malloc(sizeof(Read));
                                reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
                                reads[nbchr-1]->next->offset = lineSize;
                                reads[nbchr-1] = reads[nbchr-1]->next;
                                readNumberByChr[nbchr-1]++;
                        }
                                 			
			offset_read_in_source_file += lineSize;

			for(i=0;i<MAX_LINE_SIZE;i++){
                                currentLine[i]=0;
                        }
                        next = tokenizer(NULL, '\n', currentLine);

                        counter++;
                }

                fprintf(stderr, "rank %d ::: counter = %zu \n", rank, counter);

                for(i=0;i<nbchr;i++){
                        preadNumberByChr[0][i] += readNumberByChr[i];
                }
                free(readNumberByChr);
}




void parser_paired(char *localData, int rank, size_t start_offset, unsigned char threshold,
		int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads){

		char *currentCarac;
		char currentLine[MAX_LINE_SIZE];
		unsigned char quality;
		unsigned int i, chr, nbchr = 0, mchr;
		int lastChr = -1;
		int next;
		size_t lineSize, offset_read_in_source_file;
		size_t coord;
		size_t *readNumberByChr;
		size_t counter = 0;
		Read **reads = *preads;

		for(i=0;i<MAX_LINE_SIZE;i++){
			currentLine[i]=0;
		}

		//we take the first line *
		//before calling parsepaired, we know that localdata is at the begining of a read
		next = tokenizer(localData,'\n', currentLine);
		offset_read_in_source_file = start_offset;

		nbchr = nbchrom;
		readNumberByChr = (size_t*)calloc(nbchr, sizeof(size_t));

		while(next){

			lineSize = strlen(currentLine) + 1;

			//we update the offset in the
			//source file
			currentLine[lineSize - 1] = '\n';
			currentLine[lineSize] = '\0';

			//GO TO FLAG
			currentCarac = strstr(currentLine, "\t");
			
			//mandatory to sort by name
			*currentCarac = '\0';
			
			currentCarac++;

			//GO TO RNAME (Chr name)
			currentCarac = strstr(currentCarac+1, "\t");
			if(lastChr == (nbchr - 1))
			{
				chr = (nbchr -1);
			}
			else
			{
				chr = getChr(currentCarac, chrNames, nbchr);
			}

			//GO TO COORD
			currentCarac = strstr(currentCarac+1, "\t");


			if (parse_mode == MODE_NAME) {
				coord = strtoull(currentLine, NULL, strlen(currentLine));
				//coord = strtoull(currentLine, NULL, strlen(currentLine));
				coord = hash_name(currentLine, 16);
				//if(!rank) printf("%s => %zu\n", currentLine, coord);
				strtoull(currentCarac, &currentCarac, 10);
			}
			else {
				//TAKE COORD AND GO TO MAPQ
				coord = strtoull(currentCarac, &currentCarac, 10);
			}

			//TAKE MAPQ AND GO TO CIGAR
			quality = strtoull(currentCarac, &currentCarac, 10);

			//GO TO RNEXT
			currentCarac = strstr(currentCarac+1, "\t");
			if(currentCarac[1] == '='){
				mchr = chr;
			}
			else if(currentCarac[1] == '*'){
				mchr = (nbchr-1);
			}
			else{
				mchr = getChr(currentCarac, chrNames, nbchr);
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
			else if ((chr < (nbchr-2)) && ( mchr < (nbchr -2))){

					//we found discordant reads
					reads[nbchr-2]->next = malloc(sizeof(Read));
					reads[nbchr-2]->next->offset_source_file=offset_read_in_source_file;
					reads[nbchr-2]->next->offset = lineSize;
					reads[nbchr-2] = reads[nbchr-2]->next;
					readNumberByChr[nbchr-2]++;

			}
			else if ((chr == '*') && ( mchr < (nbchr -2))){

					//we found discordant reads with one pair unmapped
					reads[nbchr-2]->next = malloc(sizeof(Read));
					reads[nbchr-2]->next->offset_source_file=offset_read_in_source_file;
					reads[nbchr-2]->next->offset = lineSize;
					reads[nbchr-2] = reads[nbchr-2]->next;
					readNumberByChr[nbchr-2]++;
			}
			else if ((mchr == '*') && ( chr < (nbchr -2))){

					//we found discordant reads with one pair unmapped
					reads[nbchr-2]->next = malloc(sizeof(Read));
					reads[nbchr-2]->next->offset_source_file=offset_read_in_source_file;
					reads[nbchr-2]->next->offset = lineSize;
					reads[nbchr-2] = reads[nbchr-2]->next;
					readNumberByChr[nbchr-2]++;
			}
			else if ((mchr == 65535) && ( chr == 65535)){

                                        reads[nbchr-2]->next = malloc(sizeof(Read));
                                        reads[nbchr-2]->next->offset_source_file=offset_read_in_source_file;
                                        reads[nbchr-2]->next->offset = lineSize;
                                        reads[nbchr-2] = reads[nbchr-2]->next;
                                        readNumberByChr[nbchr-2]++;
                        }
                        else if ((mchr == 65535) && ( chr < (nbchr - 1))){

                                        reads[nbchr-2]->next = malloc(sizeof(Read));
                                        reads[nbchr-2]->next->offset_source_file=offset_read_in_source_file;
                                        reads[nbchr-2]->next->offset = lineSize;
                                        reads[nbchr-2] = reads[nbchr-2]->next;
                                        readNumberByChr[nbchr-2]++;
                        }
                        else if ((chr == 65535) && ( mchr < (nbchr - 1))){

                                        reads[nbchr-2]->next = malloc(sizeof(Read));
                                        reads[nbchr-2]->next->offset_source_file=offset_read_in_source_file;
                                        reads[nbchr-2]->next->offset = lineSize;
                                        reads[nbchr-2] = reads[nbchr-2]->next;
                                        readNumberByChr[nbchr-2]++;
                        }                                                       
			else{
					//we found unmapped pairs reads
					reads[nbchr-1]->next = malloc(sizeof(Read));
					reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
					reads[nbchr-1]->next->offset = lineSize;
					reads[nbchr-1] = reads[nbchr-1]->next;
					readNumberByChr[nbchr-1]++;
			}



			//we update the offset_read_in_source_file
			offset_read_in_source_file += lineSize;
			//we read the next line

			for(i=0;i<MAX_LINE_SIZE;i++){
				currentLine[i]=0;
			}
			next = tokenizer(NULL, '\n', currentLine);

			counter++;
		}

		fprintf(stderr, "rank %d ::: counter = %zu \n", rank, counter);

		for(i=0;i<nbchr;i++){
			preadNumberByChr[0][i] += readNumberByChr[i];
		}
		free(readNumberByChr);
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
