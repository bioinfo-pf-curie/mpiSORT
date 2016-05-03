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
     preWrite.c

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/


#include <stdlib.h>
#include "preWrite.h"
#include "parser.h"
#include "diffuse.h"
#include "mpi.h"

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

MPI_Offset unmappedOffset(int rank, int num_proc, size_t blockSize, size_t headerSize, int nbchr, size_t readNum){
	size_t delay = 0;
	MPI_Offset offset = headerSize;

	//arrays[nbchr-1] = saveArrays[nbchr-1]->next;
	if(rank){

		MPI_Recv(&delay, 1, MPI_LONG_LONG_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		offset = delay;
	}

	if(rank < num_proc - 1){
		delay += blockSize;

		MPI_Send(&delay, 1, MPI_LONG_LONG_INT, rank + 1, 0, MPI_COMM_WORLD);

	}

	return offset;
}

MPI_Offset discordantOffset(int rank, int num_proc, size_t blockSize, size_t headerSize, int nbchr, size_t readNum){
	size_t delay = 0;
	MPI_Offset offset = headerSize;

	//arrays[nbchr-1] = saveArrays[nbchr-1]->next;
	if(rank){

		MPI_Recv(&delay, 1, MPI_LONG_LONG_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		offset = delay;
	}

	if(rank < num_proc - 1){
		delay += blockSize;

		MPI_Send(&delay, 1, MPI_LONG_LONG_INT, rank + 1, 0, MPI_COMM_WORLD);

	}

	return offset;
}
