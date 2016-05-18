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
     diffuse.c

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/



#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "diffuse.h"
#include "merge.h"
#include "mpi.h"

void diffuse(size_t *recv_diffuse, int rank, int num_proc, char sender, size_t localReadNum, size_t* count_diffuse, size_t **send_diffuse){

	size_t j, k;
	MPI_Request *reqs[2];
	reqs[0] = 0;
	reqs[1] = 0;
	size_t readNum, rcvd = 0;
	MPI_Status stat;
	//size_t* recv_diffuse;
	//recv_diffuse = (size_t*)malloc(localReadNum*sizeof(size_t));

	if(sender){

		reqs[0] = (MPI_Request*)malloc(num_proc*sizeof(MPI_Request));
		reqs[1] = (MPI_Request*)malloc(num_proc*sizeof(MPI_Request));

		for(j = 0; j < num_proc; j++){

			if(count_diffuse[j] && j != rank){
				MPI_Isend(&count_diffuse[j], 1, MPI_LONG_LONG_INT, j, j, MPI_COMM_WORLD, &reqs[0][j]);
				MPI_Isend(send_diffuse[j], count_diffuse[j], MPI_LONG_LONG_INT, j, j, MPI_COMM_WORLD, &reqs[1][j]);
			}
		}
	}

	j = num_proc-1;

	while(rcvd < localReadNum){

		if(j == rank){

			for(k = 0; k < count_diffuse[j]; k++)
				recv_diffuse[k+rcvd] = send_diffuse[j][k];
			rcvd += count_diffuse[j];
		}

		else{

			MPI_Recv(&readNum, 1, MPI_LONG_LONG_INT, j, rank, MPI_COMM_WORLD, &stat);
			MPI_Recv(recv_diffuse+rcvd, readNum, MPI_LONG_LONG_INT, j, rank, MPI_COMM_WORLD, &stat);

			rcvd += readNum;
		}

		j = j > 3 ? j-1 : num_proc-1;
	}

	if(sender){
		for(j = 0; j < num_proc; j++){
			if(count_diffuse[j]){
				if(j != rank){

					MPI_Wait(&reqs[0][j], MPI_STATUS_IGNORE);
					MPI_Wait(&reqs[1][j], MPI_STATUS_IGNORE);
				}
				free(send_diffuse[j]);
			}
		}

		free(send_diffuse);
		free(count_diffuse);
		count_diffuse = NULL;
		free(reqs[0]);
		free(reqs[1]);
	}

	//return recv_diffuse;
}
