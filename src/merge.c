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
     merge.c

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/



#include <stdio.h>
#include "mpi.h"
#include "merge.h"
#include "format.h"
#include "diffuse.h"
#include "merge_utils.h"

char merge(int rank, int num_proc, int headerSize, size_t readNum, size_t array_max_size, size_t **count_diffusep,
		size_t ***send_diffusep, size_t *dsend[3]){

	int father, son, dest, src, last;
	int nbSons;
	size_t nbrecv;
	char sender = 0;
	size_t wtd[2]; //what to do : wtd[0] = nbrecv wtd[1] = are you the first to have nothing to do ?

	father = (rank + rank%2 - 2)/2;

	nbSons = number_sons(rank, num_proc);

	if (rank && nbSons >= 1){

		son = 2 * rank + 1;
		readNum = datarecv(dsend, son, readNum);
	}

	if (rank && nbSons == 2){

		son = 2 * rank + 2;
		readNum = datarecv(dsend, son, readNum);

	}

	if(rank > 2){

		datasend(dsend, father, readNum);

		src = rank+1 < num_proc ? rank+1 : 3;

		free(dsend[0]);
		free(dsend[1]);
		free(dsend[2]);

		do{

			MPI_Recv(wtd, 2, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			nbrecv = wtd[0];

			if(nbrecv){

				sender = 1; //if the current rank does receive data from rank 0, it has to send data during the diffusion
				*send_diffusep = recv_dispatch(nbrecv, headerSize, rank, num_proc, count_diffusep);
			}

			else if(wtd[1]){
				//receiving data wich will not be used in order to avoid any unintended blocking - see recv_dispatch()
				MPI_Recv(&headerSize, 1, MPI_LONG_LONG_INT, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}while(nbrecv);
	}

	//ranks 1 & 2 will segment the data they send to restrict the amount of data the rank 0 will have to manage at the same time
	else if(rank)
		data_pieces_send(dsend, rank, readNum, array_max_size);

	else{

		dest = num_proc-1;

		dest = data_pieces_recv(dsend, rank, num_proc, readNum);

		last = dest;

		wtd[0] = 0;
		wtd[1] = 1;

		do{
			if(dest < 4) dest = num_proc-1;
			else dest--;
			MPI_Send(wtd, 2, MPI_LONG_LONG_INT, dest, 0, MPI_COMM_WORLD);
			wtd[1] = 0;
		}while(dest != last);
	}
	return sender;
}

size_t datarecv(size_t *send[3], int son, size_t readNum){

	size_t nbrecv;
	size_t *tmp[3];
	size_t *recv[3];

	MPI_Recv(&nbrecv,1, MPI_LONG_LONG_INT,son, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	recv[0] = (size_t*)malloc(nbrecv*sizeof(size_t));
	recv[1] = (size_t*)malloc(nbrecv*sizeof(size_t));
	recv[2] = (size_t*)malloc(nbrecv*sizeof(size_t));

	MPI_Recv(recv[0],nbrecv, MPI_LONG_LONG_INT,son, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(recv[1],nbrecv, MPI_LONG_LONG_INT,son, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(recv[2],nbrecv, MPI_LONG_LONG_INT,son, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	tmp[0] = (size_t*)malloc((readNum + nbrecv) * sizeof(size_t));
	tmp[1] = (size_t*)malloc((readNum + nbrecv) * sizeof(size_t));
	tmp[2] = (size_t*)malloc((readNum + nbrecv) * sizeof(size_t));


	arrayMerge(send, recv, tmp, readNum, nbrecv);

	readNum+=nbrecv;

	free(recv[0]);
	free(recv[1]);
	free(recv[2]);
	free(send[0]);
	free(send[1]);
	free(send[2]);

	send[0] = tmp[0];
	send[1] = tmp[1];
	send[2] = tmp[2];

	return readNum;
}

void datasend(size_t *send[3], int dest, size_t readNum){

	MPI_Send(&readNum,1, MPI_LONG_LONG_INT, dest, 0, MPI_COMM_WORLD);


	MPI_Send(send[0],readNum, MPI_LONG_LONG_INT, dest, 0, MPI_COMM_WORLD);
	MPI_Send(send[1],readNum, MPI_LONG_LONG_INT, dest, 0, MPI_COMM_WORLD);
	MPI_Send(send[2],readNum, MPI_LONG_LONG_INT, dest, 0, MPI_COMM_WORLD);
}

size_t** recv_dispatch(size_t nbrecv, size_t disp, int rank, int num_proc, size_t** count_diffusep){

	size_t *recv[3], j;
	size_t readNum;
	size_t **send_diffuse;
	size_t *abs_offset;

	int dest = rank-1 > 2 ? rank-1 : num_proc-1;
	int src = rank+1 < num_proc ? rank+1 : 3;

	recv[0] = (size_t*)malloc(nbrecv*sizeof(size_t));
	recv[1] = (size_t*)malloc(nbrecv*sizeof(size_t));
	recv[2] = (size_t*)malloc(nbrecv*sizeof(size_t));
	abs_offset = (size_t*)malloc(nbrecv*sizeof(size_t));

	MPI_Recv(recv[0], nbrecv, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(recv[1], nbrecv, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(recv[2], nbrecv, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if(rank < num_proc-1)
		MPI_Recv(abs_offset, 2, MPI_LONG_LONG_INT, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	else
		abs_offset[0] = disp;

	for (j = 1; j < nbrecv; j++){
		abs_offset[j] = abs_offset[j - 1] + recv[2][j - 1];
	}

	disp = abs_offset[j - 1] + recv[2][j - 1];

	readNum = nbrecv;

	send_diffuse = format(rank, num_proc, recv, readNum, count_diffusep, abs_offset);

	if(dest != rank){
		MPI_Send(&disp, 1, MPI_LONG_LONG_INT, dest, 1, MPI_COMM_WORLD);
	}

	return send_diffuse;
}

void data_pieces_send(size_t *send[3], int rank, size_t readNum, size_t array_max_size){

	size_t rv = 1;
	size_t cmax = 0, spent = 0;
	MPI_Status stat;
	size_t s;
	char sr = 1; //is another MPI_Sendrecv needed (to avoid erroneous blocking) ?
	int j = 0;

	MPI_Send(&readNum, 1, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);

	while(spent < readNum && rv){

		j++;

		if(readNum - spent > array_max_size){
			cmax = array_max_size;
			while(cmax && send[1][spent + cmax] == send[1][spent + cmax-1])
				cmax--;
			MPI_Sendrecv(&send[1][spent + cmax], 1, MPI_LONG_LONG_INT, rank%2+1, 1, &rv, 1, MPI_LONG_LONG_INT, rank%2+1, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
		}

		else{
			cmax = readNum - spent;
			s = 0;
			MPI_Sendrecv(&s, 1, MPI_LONG_LONG_INT, rank%2+1, 1, &rv, 1, MPI_LONG_LONG_INT, rank%2+1, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

			s = 1;
			/* can't send any value to the other rank since all the data is meant to be sent to rank 0, but the tag is put to 1 in order to avoid
			 * any misunderstanding with "end of data" (can be different if rv < send[1][readNum-1]
			 */
		}

		if(rv){
			if(send[1][spent] > rv) //if reads are grouped, will greatly shorten the amount of comparisons
				s = 0; //no data will be sent to rank 0 in this iteration

			else{
				while(cmax && send[1][spent + cmax-1] >= rv){
					cmax--;
				}
			}
		}

		else if(stat.MPI_TAG)
			rv = 1; // since the other rank is supposed to send all of its data to rank 0, rv isn't useful. However, there is data to send

		else
			sr = 0;

		if(s){
			MPI_Send(&cmax, 1, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);

			MPI_Send(send[0] + spent, cmax, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Send(send[1] + spent, cmax, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Send(send[2] + spent, cmax, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);

			spent += cmax;
		}

		else{
			MPI_Send(&s, 1, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);
		}
	}

	if(sr){

		s = 0;

		MPI_Sendrecv(&s, 1, MPI_LONG_LONG_INT, rank%2+1, 0, &rv, 1, MPI_LONG_LONG_INT, rank%2+1, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

		//tag = 0 and value (s) is 0 : end of data
	}

	else while(spent < readNum){

		j++;

		if(readNum - spent > array_max_size){
			cmax = array_max_size;
			while(cmax && send[1][spent + cmax] == send[1][spent + cmax-1])
				cmax--;
		}

		else
			cmax = readNum - spent;

		MPI_Send(&cmax, 1, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);

		MPI_Send(send[0] + spent, cmax, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(send[1] + spent, cmax, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(send[2] + spent, cmax, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);

		spent += cmax;
	}

	free(send[0]);
	free(send[1]);
	free(send[2]);
}

int data_pieces_recv(size_t *send[3], int rank, int num_proc, size_t readNum){

	size_t *recv1[3], *recv2[3], *tmp[3];
	size_t spent = 0, last_spent = 0;
	size_t readNum1, readNum2, tmp1, tmp2;
	int k = num_proc-1;
	tmp1 = tmp2 = 1;
	int j = 0;

	MPI_Recv(&readNum1, 1, MPI_LONG_LONG_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&readNum2, 1, MPI_LONG_LONG_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	while(readNum1 && readNum2){

		j++;

		if(tmp1)MPI_Recv(&tmp1, 1, MPI_LONG_LONG_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if(tmp2)MPI_Recv(&tmp2, 1, MPI_LONG_LONG_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		if(tmp1){

			recv1[0] = (size_t*)malloc(tmp1*sizeof(size_t));
			recv1[1] = (size_t*)malloc(tmp1*sizeof(size_t));
			recv1[2] = (size_t*)malloc(tmp1*sizeof(size_t));

			MPI_Recv(recv1[0], tmp1, MPI_LONG_LONG_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(recv1[1], tmp1, MPI_LONG_LONG_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(recv1[2], tmp1, MPI_LONG_LONG_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		if(tmp2){

			recv2[0] = (size_t*)malloc(tmp2*sizeof(size_t));
			recv2[1] = (size_t*)malloc(tmp2*sizeof(size_t));
			recv2[2] = (size_t*)malloc(tmp2*sizeof(size_t));

			MPI_Recv(recv2[0], tmp2, MPI_LONG_LONG_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(recv2[1], tmp2, MPI_LONG_LONG_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(recv2[2], tmp2, MPI_LONG_LONG_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		tmp[0] = (size_t*)malloc((tmp1+tmp2+readNum-spent)*sizeof(size_t));
		tmp[1] = (size_t*)malloc((tmp1+tmp2+readNum-spent)*sizeof(size_t));
		tmp[2] = (size_t*)malloc((tmp1+tmp2+readNum-spent)*sizeof(size_t));

		readNum1 -= tmp1;
		readNum2 -= tmp2;

		spent = bucket_arrayMerge(recv1, recv2, send, tmp, tmp1, tmp2, readNum, (readNum1 || readNum2), last_spent);

		free(recv1[0]);
		free(recv1[1]);
		free(recv1[2]);
		free(recv2[0]);
		free(recv2[1]);
		free(recv2[2]);

		datasend(tmp, k, tmp1+tmp2+spent-last_spent);

		last_spent = spent;
		k = k > 3 ? k-1 : num_proc-1;

		free(tmp[0]);
		free(tmp[1]);
		free(tmp[2]);
	}

	while(readNum1){

		MPI_Recv(&tmp1, 1, MPI_LONG_LONG_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		if(tmp1){

			recv1[0] = (size_t*)malloc(tmp1*sizeof(size_t));
			recv1[1] = (size_t*)malloc(tmp1*sizeof(size_t));
			recv1[2] = (size_t*)malloc(tmp1*sizeof(size_t));

			MPI_Recv(recv1[0], tmp1, MPI_LONG_LONG_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(recv1[1], tmp1, MPI_LONG_LONG_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(recv1[2], tmp1, MPI_LONG_LONG_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}


		tmp[0] = (size_t*)malloc((tmp1+readNum)*sizeof(size_t));
		tmp[1] = (size_t*)malloc((tmp1+readNum)*sizeof(size_t));
		tmp[2] = (size_t*)malloc((tmp1+readNum)*sizeof(size_t));

		readNum1 -= tmp1;

		spent = bucket_arrayMerge(recv1, NULL, send, tmp, tmp1, 0, readNum, readNum1, last_spent);

		free(recv1[0]);
		free(recv1[1]);
		free(recv1[2]);

		readNum -= spent + last_spent;

		datasend(tmp, k, tmp1+spent-last_spent);

		last_spent = spent;
		k = k > 3 ? k-1 : num_proc-1;

		free(tmp[0]);
		free(tmp[1]);
		free(tmp[2]);
	}

	while(readNum2){

		MPI_Recv(&tmp2, 1, MPI_LONG_LONG_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		if(tmp2){

			recv2[0] = (size_t*)malloc(tmp2*sizeof(size_t));
			recv2[1] = (size_t*)malloc(tmp2*sizeof(size_t));
			recv2[2] = (size_t*)malloc(tmp2*sizeof(size_t));

			MPI_Recv(recv2[0], tmp2, MPI_LONG_LONG_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(recv2[1], tmp2, MPI_LONG_LONG_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(recv2[2], tmp2, MPI_LONG_LONG_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		tmp[0] = (size_t*)malloc((tmp2+readNum)*sizeof(size_t));
		tmp[1] = (size_t*)malloc((tmp2+readNum)*sizeof(size_t));
		tmp[2] = (size_t*)malloc((tmp2+readNum)*sizeof(size_t));

		readNum2 -= tmp2;

		spent = bucket_arrayMerge(NULL, recv2, send, tmp, 0, tmp2, readNum, readNum2, last_spent);

		free(recv2[0]);
		free(recv2[1]);
		free(recv2[2]);

		readNum -= spent + last_spent;

		datasend(tmp, k, tmp2+spent-last_spent);
		last_spent = spent;

		k = k > 3 ? k-1 : num_proc-1;

		free(tmp[0]);
		free(tmp[1]);
		free(tmp[2]);
	}

	free(send[0]);
	free(send[1]);
	free(send[2]);

	k = k < num_proc-1 ? k+1 : 3;

	return k;
}
