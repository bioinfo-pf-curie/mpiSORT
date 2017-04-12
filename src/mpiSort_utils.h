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
     mpiSort_utils.h

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void get_coordinates_and_offset_source_and_size_and_free_reads(
		int rank,
		int *local_read_rank,
		size_t *coordinates,
		size_t* offset,
		int* size,
		Read* data_chr,
		int local_readNum
		);

size_t init_coordinates_and_size(
		int rank,
		int *local_reads_rank,
		size_t *local_reads_index,
		size_t* coordinates,
		int* size,
		Read* data_chr,
		int local_readNum
		);


void chosen_split_rank_gather_size_t(
		MPI_Comm split_comm,
		int rank,
		int num_proc,
		int master,
		size_t size,
		size_t *size_per_jobs,
		size_t *start_size_per_job,
		size_t *all_data,
		size_t *data,
		size_t start_index
		);
