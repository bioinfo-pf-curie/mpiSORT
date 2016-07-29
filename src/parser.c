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
     parser.c

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/



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

unsigned int find_header(char *localData, int rank, size_t *unmappedSize, int *pnbchr, char **pheader, char ***pchrNames){
	char *currentCarac;
	char *header;
	char **chrNames = NULL;
	char currentLine[MAX_LINE_SIZE];
	unsigned int i, nbchr = 0;
	int next;
	size_t headerSize, lineSize;
	size_t j;

	//we take the first line
	next = tokenizer(localData,'\n', currentLine);

	headerSize = 0; //offset size of header to replace offset of proc 0
	lineSize = 0;//use to update headerSize

	chrNames = (char**)malloc(sizeof(char*));

	char *str = "@HD\tVN:1.0\tSO:coordinate";

	//RANK 0 WILL CALCULATE SIZE HEADER AND FIND CHRNAME
	if(rank == 0){

		j = 0;
		header = (char*)malloc(sizeof(char));
		lineSize = strlen(str) + 1;
		headerSize += lineSize;
		header = realloc(header, (headerSize+1)*sizeof(char));
		while(*str){
				header[j++] = *str++;
		}
		header[j++] = '\n';

		//we look headers lines
		//we add header string into variable 'header' and add find chrNames
		while (next && currentLine[0] == '@'){

			currentCarac = currentLine;

			lineSize = strlen(currentLine) + 1;
			headerSize += lineSize;

			// in previous version we use realloc
			header = realloc(header, (headerSize+1)*sizeof(char));

			//add currentline into header
			while(*currentCarac){
				header[j++] = *currentCarac++;
			}

			header[j++] = '\n';
			header[j] = 0;
			*pheader = header;

			//find chromosome name
			if(currentLine[1] == 'S' && currentLine[2] == 'Q'){

				currentCarac = strtok(currentLine, "\t");
				currentCarac = strtok(NULL, "\t");
				currentCarac += 3;

				chrNames[nbchr++] = strdup(currentCarac);

				chrNames = (char**)realloc(chrNames, (nbchr+1)*sizeof(char*));

			}

			next = tokenizer(NULL, '\n', currentLine);
		}

		//header[j++] = 0;

		//we add UNMAPPED to chrNames
		chrNames[nbchr++] = strdup(UNMAPPED);

		//we add DISCORDANT to chrNames
		chrNames[nbchr++] = strdup(DISCORDANT);

	}

	/*************************
	 * GIVE THE HEADER SIZE TO ALL
	 *************************/
	if(rank==0){
		fprintf(stderr, "HEADERSIZE : %zu\n",headerSize);
		fprintf(stderr, "%s\n",header);
		fprintf(stderr, "NB CHR : %d\n",nbchr);
	}
	MPI_Bcast(&headerSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	*unmappedSize = headerSize;


	/**********************
	 * GIVE CHR NAMES TO ALL
	 *************************/

	//1 - every process will know nbchr
	MPI_Bcast(&nbchr, 1, MPI_INT, 0,MPI_COMM_WORLD);
	*pnbchr = nbchr;

	if(nbchr!=0){
		size_t size_chrName[nbchr];//list of size of chromosome name
		size_t size_all_chrname = 0;
		for(i=0;i<nbchr;i++){
			size_chrName[i]=0;
		}

		//2 - rank 0 will set the size of each chrname
		if(rank ==0){
			for(i=0;i<nbchr;i++){
				size_chrName[i]=strlen(chrNames[i]);
				size_all_chrname+=size_chrName[i];
			}
			size_all_chrname++;
		}

		//3 - Other process will prepare buff for each chromosome name
		MPI_Bcast(size_chrName,nbchr,MPI_UNSIGNED_LONG,0,MPI_COMM_WORLD);
		if(rank){
			chrNames = (char**)malloc((nbchr)*sizeof(char*));

			for(i=0;i<nbchr;i++){
				chrNames[i]=(char*)malloc(size_chrName[i]*sizeof(char)+1);
				chrNames[i][size_chrName[i]]='\0';

				*(chrNames[i]) = 0;
				size_all_chrname += size_chrName[i];
			}
			size_all_chrname++;
		}

		//4 - We will send all read names in one buffer
		char * buff_chrNames =(char*)malloc(size_all_chrname*sizeof(char)+1);

		for(i=0;i<size_all_chrname;i++){
			buff_chrNames[i]=0;
		}

		buff_chrNames[size_all_chrname]=0;

		if(rank==0){
			for(i=0;i<nbchr;i++){
				strcat(buff_chrNames,chrNames[i]);
			}
		}

		//5 - Other rank will add chr names
		MPI_Bcast(buff_chrNames,(int)size_all_chrname,MPI_CHAR,0,MPI_COMM_WORLD);
		if(rank){
			int offset_chrname = 0;
			for(i=0;i<nbchr;i++){
				int n = 0;
				for(n=0;n<size_chrName[i];n++){
					chrNames[i][n]=buff_chrNames[offset_chrname+n];
				}
				chrNames[i][size_chrName[i]]='\0';
				offset_chrname+=size_chrName[i];
			}
		}

		free(buff_chrNames);
		*pchrNames = chrNames;
	}
	return(headerSize);
}

size_t * init_goff(MPI_File mpi_filed,unsigned int headerSize,size_t fsize,int numproc,int rank){
	size_t * goff =(size_t*)calloc((size_t)(numproc+1), sizeof(size_t));
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
		MPI_File_read_at(mpi_filed, (MPI_Offset)goff[i], current_line, 1000, MPI_CHAR, &status);
		j=0;
		while(j<fsize && current_line[j] != '\n'){
			j++;
		}
		goff[i]+=(j+1);
		free(current_line);
	}

	return goff;
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
				if(!rank) printf("%s => %zu\n", currentLine, coord);
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
					if(quality > threshold){

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
					reads[nbchr-1]->next = malloc(sizeof(Read));
					reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
					reads[nbchr-1]->next->offset = lineSize;
					reads[nbchr-1] = reads[nbchr-1]->next;
					readNumberByChr[nbchr-1]++;

			}
			else if ((chr == '*') && ( mchr < (nbchr -2))){

					//we found discordant reads with one pair unmapped
					reads[nbchr-1]->next = malloc(sizeof(Read));
					reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
					reads[nbchr-1]->next->offset = lineSize;
					reads[nbchr-1] = reads[nbchr-1]->next;
					readNumberByChr[nbchr-1]++;
			}
			else if ((mchr == '*') && ( chr < (nbchr -2))){

					//we found discordant reads with one pair unmapped
					reads[nbchr-1]->next = malloc(sizeof(Read));
					reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
					reads[nbchr-1]->next->offset = lineSize;
					reads[nbchr-1] = reads[nbchr-1]->next;
					readNumberByChr[nbchr-1]++;
			}

			else{
					//we found unmapped pairs reads
					reads[nbchr-2]->next = malloc(sizeof(Read));
					reads[nbchr-2]->next->offset_source_file=offset_read_in_source_file;
					reads[nbchr-2]->next->offset = lineSize;
					reads[nbchr-2] = reads[nbchr-2]->next;
					readNumberByChr[nbchr-2]++;
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
		//free(readNumberByChr);
}

void parser_single(char *localData, int rank, size_t start_offset, unsigned char threshold,
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
			//TAKE COORD AND GO TO MAPQ
			coord = strtoull(currentCarac, &currentCarac, 10);

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
					if(quality > threshold){

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
					reads[nbchr-1]->next = malloc(sizeof(Read));
					reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
					reads[nbchr-1]->next->offset = lineSize;
					reads[nbchr-1] = reads[nbchr-1]->next;
					readNumberByChr[nbchr-1]++;

			}
			else if ((chr == '*') && ( mchr < (nbchr -2))){

					//we found discordant reads with one pair unmapped
					reads[nbchr-1]->next = malloc(sizeof(Read));
					reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
					reads[nbchr-1]->next->offset = lineSize;
					reads[nbchr-1] = reads[nbchr-1]->next;
					readNumberByChr[nbchr-1]++;
			}
			else if ((mchr == '*') && ( chr < (nbchr -2))){

					//we found discordant reads with one pair unmapped
					reads[nbchr-1]->next = malloc(sizeof(Read));
					reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
					reads[nbchr-1]->next->offset = lineSize;
					reads[nbchr-1] = reads[nbchr-1]->next;
					readNumberByChr[nbchr-1]++;
			}

			else{
					//we found unmapped pairs reads
					reads[nbchr-2]->next = malloc(sizeof(Read));
					reads[nbchr-2]->next->offset_source_file=offset_read_in_source_file;
					reads[nbchr-2]->next->offset = lineSize;
					reads[nbchr-2] = reads[nbchr-2]->next;
					readNumberByChr[nbchr-2]++;
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
	return i-1;
}
