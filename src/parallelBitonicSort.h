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
     parallelBitonicSort.h

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/

#include <mpi.h>

//typedef size_t size_t;
void Generate_local_list(
		size_t list_size,
		size_t local_list[],
		size_t local_index[]
		);
void Local_sort(
		size_t list_size,
		size_t local_keys[],
		size_t local_index[]
		);
int Key_compare(
		const size_t* p,
		const size_t* q
		);
void Par_bitonic_sort_incr(
		size_t list_size,
		size_t local_list[],
		size_t index_list[],
        int proc_set_size,
        int rank
        );
void Par_bitonic_sort_decr(
		size_t list_size,
		size_t local_list[],
		size_t index_list[],
        int proc_set_size,
        int rank
        );
void Merge_split(
		size_t list_size,
		size_t local_list[],
		size_t local_index[],
		int which_keys,
		int partner,
		int my_rank
		);
void Merge_list_low(
		size_t  list_size,
		size_t  list_key[],
		size_t  list_index[],
		size_t  list_tmp_key[],
		size_t  list_tmp_index[]
		);
void Merge_list_high(
		size_t  list_size,
		size_t  list_key[],
		size_t  list_index[],
		size_t  list_tmp_key[],
		size_t  list_tmp_index[]
		);
int bitonic_qksort(void *data,
		size_t size,
		size_t esize,
		size_t i,
		size_t k,
		int (*compare)(const void *key1, const void *key2)
		);

int bitonic_partition(void *data,
		size_t esize,
		size_t i,
		size_t k,
		int (*compare)(const void *key1, const void *key2)
		);

void ParallelBitonicSort(
		MPI_Comm split_comm,
		int my_rank,
		int dimension,
		size_t *local_list,
		size_t *local_index,
		size_t list_size,
		size_t zero_padding
		);
