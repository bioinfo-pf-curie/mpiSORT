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
     parabitonicsort.h

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/



#include <mpi.h>
//typedef size_t size_t;
void Generate_local_list(int list_size, size_t local_list[], size_t local_index[]);
void Local_sort(size_t list_size, size_t local_keys[], size_t local_index[]);
int Key_compare(const size_t* p, const size_t* q);
int log_base2(int x);
void Par_bitonic_sort_incr(int list_size, size_t local_list[], size_t index_list[],
        int proc_set_size);
void Par_bitonic_sort_decr(int list_size, size_t local_list[], size_t index_list[],
        int proc_set_size);
void Merge_split(size_t list_size, size_t local_list[], size_t local_index[],int which_keys, int partner, int my_rank);
void Merge_list_low(size_t list_size, size_t  list_key[], size_t  list_index[], size_t list_tmp_key[], size_t  list_tmp_index[]);
void Merge_list_high(size_t list_size, size_t  list_key[], size_t  list_index[], size_t list_tmp_key[], size_t  list_tmp_index[]);
int bitonic_qksort(void *data, size_t size, size_t esize, size_t i, size_t k, int (*compare)(const void *key1, const void *key2));
int bitonic_partition(void *data, size_t esize, size_t i, size_t k, int (*compare)(const void *key1, const void *key2));
void ParallelBitonicSort(MPI_Comm split_comm, int my_rank, int dimension, size_t *local_list, size_t *local_index, size_t list_size, size_t zero_padding);
