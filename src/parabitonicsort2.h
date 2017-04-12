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
     parabitonicsort.h

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
void Generate_local_list2(
		size_t list_size,
		size_t local_list[]
		);
void Local_sort2(
		size_t list_size,
		size_t local_keys[],
		int    local_keys1[],
		int    local_keys2[],
		size_t local_keys3[],
		int local_keys4[]
		);

int Key_compare2(const size_t* p, const size_t* q);
void Par_bitonic_sort_incr2(
		size_t list_size,
		size_t local_list[],
		int    local_list1[],
		int    local_list2[],
		size_t local_list3[],
		int local_list4[],
		int    proc_set_size,
        int    rank
        );
void Par_bitonic_sort_decr2(
		size_t list_size,
		size_t local_list[],
		int    local_list1[],
		int    local_list2[],
		size_t local_list3[],
		int local_list4[],
		int    proc_set_size,
        int    rank
        );
void Merge_split2(
		size_t list_size,
		size_t local_list[],
		int    local_list1[],
		int    local_list2[],
		size_t local_list3[],
		int local_list4[],
		int    which_keys,
		int    partner,
		int    my_rank
		);
void Merge_list_low2(
		size_t  list_size,
		size_t  list_key[],
		int     list_key1[],
		int     list_key2[],
		size_t  list_key3[],
		int  list_key4[],
		size_t  list_tmp_key[],
		int     list_tmp_key1[],
		int     list_tmp_key2[],
		size_t  list_tmp_key3[],
		int  list_tmp_key4[]
		);
void Merge_list_high2(
		size_t  list_size,
		size_t  list_key[],
		int     list_key1[],
		int     list_key2[],
		size_t  list_key3[],
		int     list_key4[],
		size_t  list_tmp_key[],
		int     list_tmp_key1[],
		int     list_tmp_key2[],
		size_t  list_tmp_key3[],
		int     list_tmp_key4[]
		);
int bitonic_qksort2(
		void   *data,
		size_t size,
		size_t esize,
		size_t i,
		size_t k,
		int    (*compare)(const void *key1, const void *key2)
		);
int bitonic_partition2(
		void   *data,
		size_t esize,
		size_t i,
		size_t k,
		int    (*compare)(const void *key1, const void *key2)
		);
void ParallelBitonicSort2(
		MPI_Comm split_comm,
		int      my_rank,
		int      dimension,
		size_t   *local_list,
		int      *local_list1,
		int      *local_list2,
		size_t   *local_list3,
		int      *local_list4,
		size_t   list_size
		);
