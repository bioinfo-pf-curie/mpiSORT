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
     format.c

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/



#include <stdio.h>
#include <stdlib.h>
#include "merge.h"
#include "parser.h"
#include "diffuse.h"
#include "format.h"

size_t** format(int rank, int num_proc, size_t* recv[3], size_t readNum, size_t** count_diffusep, size_t *abs_offset){

	/*FIXME: For now, if there are more than 4 processes, one process can't format more than once.
	 * 		 In order to change that, we need to add one dimension to (send-count-iterator)_diffuse, and the first dimension will be used
	 *		 for the number of times the same process will format, even though this is quite unlikely.
	 * */

	size_t j;
	size_t k;
	size_t rankDest;
	size_t *iterator_diffuse = NULL, *count_diffuse = NULL;
	size_t **send_diffuse = NULL;

	if(!*count_diffusep || num_proc > 4){

		count_diffuse = (size_t*)calloc(num_proc, sizeof(size_t));
		iterator_diffuse = (size_t*)calloc(num_proc, sizeof(size_t));

		countNumRead(num_proc, recv, readNum, count_diffuse);

		send_diffuse = (size_t**)malloc(num_proc*sizeof(size_t*));

		if(send_diffuse != NULL){

			for(j = 0; j < num_proc; j++){

				if(count_diffuse[j] != 0)
					send_diffuse[j] = (size_t*)malloc(count_diffuse[j]*sizeof(size_t));
			}

			for(j = 0; j < readNum; j++){

				rankDest = recv[0][j];
				k = iterator_diffuse[rankDest];

				send_diffuse[rankDest][k] = abs_offset[j];

				iterator_diffuse[rankDest]++;
			}
		}
	}

	else{
		countNumRead(num_proc, recv, readNum, count_diffuse);

		for(j = 0; j < num_proc; j++){

			if(count_diffuse[j] != 0)
				send_diffuse[j] = (size_t*)realloc(send_diffuse[j], count_diffuse[j]*sizeof(size_t));
		}

		for(j = 0; j < readNum; j++){

			rankDest = recv[0][j];
			k = iterator_diffuse[rankDest];

			send_diffuse[rankDest][k] = abs_offset[j];

			iterator_diffuse[rankDest]++;
		}
	}

	free(abs_offset);
	free(recv[0]);
	free(recv[1]);
	free(recv[2]);
	free(iterator_diffuse);
	*count_diffusep = count_diffuse;

	return send_diffuse;
}

/**
 * \fn void countNumRead(int num_proc)
 * \brief count the number of reads of each process
 * \param num_proc is the processes total number
 */

void countNumRead(int num_proc, size_t* recv[3], size_t readNum, size_t *count_diffuse){

	size_t j;
	size_t rankDest;

	for(j = 0; j < readNum; j++){
		rankDest = recv[0][j];
		count_diffuse[rankDest]++;
	}
}
