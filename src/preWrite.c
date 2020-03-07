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
     preWrite.c

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

#include <stdlib.h>

#include <mpi.h>

#include "preWrite.h"
#include "writeUtils.h"
#include "parser.h"

void preWrite(Read* chr, size_t readNum, char** data, int* size){
	size_t j;
	size_t k=0;
	size_t dataSize=0;
	char* s;
	void* tmp;

	for(j = 0; j < readNum; j++){

		size[j] = strlen(chr->string);
		dataSize += size[j];
		s = chr->string;

		*data = (char*)realloc(data, (dataSize+1)*sizeof(char));

		while(*s){
			(*data)[k++] = *s++;
		}

		free(chr->string);
		tmp = chr;
		chr = chr->next;
		free(tmp);
	}
	(*data)[k] = 0;
}

MPI_Offset startOffset(int rank, int num_proc, size_t blockSize, size_t headerSize, int nbchr, size_t readNum, MPI_Comm comm){
	size_t delay = 0;
	MPI_Offset offset = headerSize;

	//arrays[nbchr-1] = saveArrays[nbchr-1]->next;
	if(rank){

		MPI_Recv(&delay, 1, MPI_LONG_LONG_INT, rank - 1, 0, comm, MPI_STATUS_IGNORE);

		offset = delay;
	}

	if(rank < num_proc - 1){
		delay += blockSize;

		MPI_Send(&delay, 1, MPI_LONG_LONG_INT, rank + 1, 0, comm);

	}

	return offset;
}
