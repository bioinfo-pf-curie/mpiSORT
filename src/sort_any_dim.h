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
    sort_any_dim.h

   	Authors:

    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <mpi.h>

#include "qksort.h"
#include "write_utils.h"
#include "parabitonicsort.h"


void parallel_sort_any_dim(						//dimensions for parabitonic
		int dimensions,
		size_t local_readNum,
		int split_rank,
		int split_size,
		Read **reads,
		int i, 									//chromosom number
		int chosen_split_rank,
		MPI_Comm split_comm,
		size_t *localReadNumberByChr,
		char *local_data,
		char *file_name,
		char *output_dir,
		MPI_Info finfo,
		int compression_level,
		size_t total_reads_by_chr,
		size_t start_offset_in_file,
		size_t headerSize,
		char* header,
		char *chrName,
		MPI_File mpi_file_split_comm
		);
