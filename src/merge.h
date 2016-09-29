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
#include "parser.h"

char merge(int rank, int num_proc, int headerSize, size_t readNum, size_t array_max_size,
		size_t **count_diffusep, size_t ***send_diffusep, size_t *send[3]);
size_t datarecv(size_t *send[3], int src, size_t readNum);
void datasend(size_t *send[3], int dest, size_t readNum);
size_t** recv_dispatch(size_t nbrecv, size_t disp, int rank, int num_proc, size_t **count_diffusep);
void data_pieces_send(size_t *send[3], int rank, size_t readNum, size_t array_max_size);
int data_pieces_recv(size_t *send[3], int rank, int num_proc, size_t readNum);
