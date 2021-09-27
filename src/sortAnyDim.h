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
    sortAnyDim.h

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

#include "qkSort.h"
#include "writeUtils.h"
#include "parallelBitonicSort.h"


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
		int mpi_file_split_comm,
		int uniq_chr
		);
