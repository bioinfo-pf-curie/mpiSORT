/*
   mpiSORT
   Copyright (C) 2016-2017 Institut Curie / Institut Pasteur

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
#include "write_utils.h"
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
