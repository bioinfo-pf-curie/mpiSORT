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
     bufferized_read.c

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/


#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include "bufferized_read.h"
#include <assert.h>

#define DEFAULT_INBUF_SIZE (512*1024*1024)
//#define DEFAULT_INBUF_SIZE (512)
#define DEFAULT_OUTBUF_SIZE (64*1024)

void SafeFree(void **pp){

	if (pp !=NULL && *pp!=NULL){
		free(*pp);
		*pp=NULL;
	}
}

size_t StrLength(char *str){

	size_t length=0;
	while (*(str++)){
		length++;
	}
	return length;
}

void read_mpi2(MPI_File f, MPI_Comm comm, size_t start,size_t locals, char **final_data, int rank) {

		/*
		 * first we have the case th locals to read
		 * fit the MPI_File read buffer.
		 * The second case is when locals does not fit the MPI_File read buffer.
		 */

		//the final data will recieve all the data read
		*final_data = (char *)malloc(locals + 1);
		//char* interm_buffer = (char *)malloc((DEFAULT_INBUF_SIZE + 1));
		//assert(interm_buffer != NULL);
		char* interm_buffer2 = (char *)malloc((DEFAULT_INBUF_SIZE + 1));
		assert(interm_buffer2 != NULL);
		*(*final_data + locals)='\0';

		//interm_buffer[DEFAULT_INBUF_SIZE]='\0';
		interm_buffer2[DEFAULT_INBUF_SIZE]='\0';

		if (final_data == 0){
			fprintf(stderr,"Proc rank %d in ks_init_mpi problem during malloc of final_data \n", rank);
		}
		/*
		if (interm_buffer == 0){
			fprintf(stderr,"Proc rank %d in ks_init_mpi problem during malloc of interm_buffer \n", rank);
		}
		*/
		if (interm_buffer2 == 0){
			fprintf(stderr,"Proc rank %d in ks_init_mpi problem during malloc of interm_buffer2 \n", rank);
		}


		/*
		 * we create a intermediate buffer to
		 * because MPI_read_at_all cannot hold the
		 * total amount of data
		 */

		fprintf(stderr,"Proc rank %d locals =  %zu \n", rank, locals);
		fprintf(stderr,"Proc rank %d start =  %zu \n", rank, start);
		MPI_Status status;
		size_t pos_in_final_data = 0;
		size_t start1 = start;

		size_t next_start = start + DEFAULT_INBUF_SIZE;
		size_t u0;

		while( next_start < (start + locals) ){

			//with collective read cannot open 2 separate file ????
			//MPI_File_set_view(f, start1, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
			MPI_File_read_at_all(f, (MPI_Offset)start1, interm_buffer2, DEFAULT_INBUF_SIZE, MPI_CHAR, &status);

			for (u0 = 0; u0 < DEFAULT_INBUF_SIZE; u0++){
				*(*final_data + pos_in_final_data)= interm_buffer2[u0];
				pos_in_final_data +=1;
			}

			//TODO : shall we realloc???
			//interm_buffer2 = realloc(interm_buffer2, DEFAULT_INBUF_SIZE + 1);

			if (interm_buffer2 == 0){
				fprintf(stderr,"Proc rank %d !!!!!!!in ks_init_mpi problem of realloc !!!!!!!\n", rank);
			}

			start1 += DEFAULT_INBUF_SIZE;
			next_start += DEFAULT_INBUF_SIZE;

		}

		free(interm_buffer2);

		//we read the remaining data
		//we go backward because of
		//the last increment

		if (((start + locals) - start1) > 0){

			size_t remain_data_sz = (start + locals) - start1;

			char *interm_buffer3 = (char *)malloc((remain_data_sz + 1));
			interm_buffer3[remain_data_sz] = '\0';

			//with collective read cannot open 2 separate file ????
			//MPI_File_set_view(f, start1, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

			MPI_File_read_at_all(f, (MPI_Offset)start1, interm_buffer3, (MPI_Offset)remain_data_sz, MPI_CHAR, &status);

			// we copy interm buffer
			// in final_buffer
			u0 = 0;
			while ((u0 < remain_data_sz) && (pos_in_final_data <= locals)){

				*(*final_data + pos_in_final_data) = interm_buffer3[u0];
				pos_in_final_data +=1;
				u0++;
			}
			free(interm_buffer3);
		}

		fprintf(stderr,"Proc rank %d finish in read_mpi2  \n", rank);

}


