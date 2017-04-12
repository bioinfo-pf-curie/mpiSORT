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
     parabitonicsort.c

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* parallel_bitonic.c -- parallel bitonic sort of randomly generated list
 *     of integers
 *
 * Input:
 *     n: the global length of the list -- must be a power of 2.
 *
 * Output:
 *     The sorted list.
 *
 * Notes:
 *     1.  Assumes the number of processes p = 2^d and p divides n.
 *     2.  The lists are statically allocated -- size specified in MAX.
 *     3.  Keys are in the range 0 -- KEY_MAX-1.
 *     4.  Implementation can be made much more efficient by using
 *         pointers and avoiding re-copying lists in merges.
 *
 * See Chap 14, pp. 320 & ff. in PPMPI.
 *
 * https://www.eecis.udel.edu/~saunders/courses/372/01f/ppmpi_c/chap08/
 * https://www.eecis.udel.edu/~saunders/courses/372/01f/ppmpi_c/chap14a/
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "compat.h"

#include "parabitonicsort.h"
#include "parabitonicsort2.h"

// we limit to 1gb per proc

#define LOW 0
#define HIGH 1

// the number of key max is 1gb
#define key_mpi_t MPI_LONG_LONG_INT
#define index_mpi_t MPI_LONG_LONG_INT

static MPI_Comm COMM_WORLD;

/********************************************************************/

size_t *base_arr2;


static int compare_size_t_V2(const void *a, const void *b){

	if (*(const size_t *)a > *(const size_t *)b)
		return 1;
	else if (*(const size_t *)a < *(const size_t *)b)
		return -1;
	else
		return 0;

}

static int compare_size_t(const void *a, const void *b){

	 size_t aa = *(size_t *)a, bb = *(size_t *)b;
	 //fprintf(stderr, "rank %d :::::[MPIBITONIC2] we call parallel bitonic sort \n", split_rank);

	 if (base_arr2[aa] > base_arr2[bb])
		return 1;
	else if (base_arr2[aa] < base_arr2[bb])
		return -1;
	else
		return 0;
}


void ParallelBitonicSort(
		MPI_Comm split_comm,
		int my_rank,
		int dimension,
		size_t *local_list,		//stand for coordinates
		size_t *local_index, 	//stand for index
		size_t list_size,
		size_t zero_padding

) {

	// In this version of bitonic there is no more index.
	// Everything is sorted according
	// the local_list vector.

    int       proc_set_size;
    unsigned  and_bit;
    size_t k = 0;
    COMM_WORLD = split_comm;

    if (my_rank < (dimension - 1)){
    	for (k = 0; k < zero_padding; k++){
    		local_list[list_size - k - 1] = 0;
    	}
    }

    Local_sort(list_size, local_list, local_index);

    //we check the local_list is sorted
    for (k = 0; k < (list_size - 1); k++){
    	assert(local_list[k] <= local_list[k + 1]);
    }

    /* and_bit is a bitmask that, when "anded" with  */
    /* my_rank, tells us whether we're working on an */
    /* increasing or decreasing list                 */
    for (proc_set_size = 2, and_bit = 2; proc_set_size <= dimension;
    		proc_set_size = proc_set_size*2, and_bit = and_bit << 1){

        if ((my_rank & and_bit) == 0){

            Par_bitonic_sort_incr(
            		list_size,
            		local_list,
            		local_index,
            		proc_set_size,
            		my_rank
            );
        }
        else{

            Par_bitonic_sort_decr(
            		list_size,
            		local_list,
            		local_index,
            		proc_set_size,
            		my_rank
            );
        }
    }

    for (k = 0; k < (list_size - 1); k++){
      	assert(local_list[k] <= local_list[k + 1]);
    }
}


/*********************************************************************/
void Local_sort(
         size_t   list_size     /* in     */,
         size_t  *local_keys    /* in/out */,
         size_t  *local_keys1   /* in/out */
         ) {

	//we create an index vector
	size_t *local_keys_temp   	 = malloc(sizeof(size_t)*list_size);
	size_t *local_keys_temp1     = malloc(sizeof(size_t)*list_size);

	assert(local_keys_temp);
	assert(local_keys_temp1);

	size_t *index_vector = (size_t *)malloc(sizeof(size_t)*list_size);

	size_t j = 0;

	for(j = 0; j < list_size; j++){
			index_vector[j] = j;
	}

	base_arr2 = local_keys;
	bitonic_qksort(index_vector, list_size, sizeof(size_t), 0, list_size - 1, compare_size_t);

	//then we apply loac index to local_keys
	for(j = 0; j < list_size; j++){
		local_keys_temp[j]  = local_keys[index_vector[j]];
		local_keys_temp1[j] = local_keys1[index_vector[j]];
	}

	for(j = 0; j < list_size; j++){
		local_keys[j]  = local_keys_temp[j];
		local_keys1[j] = local_keys_temp1[j];
	}

	free(index_vector);
	free(local_keys_temp);
	free(local_keys_temp1);
	malloc_trim(0);
}


/*********************************************************************/
int Key_compare(const size_t* p, const size_t* q) {

    if (*p < *q)
        return -1;
    else if (*p == *q)
        return 0;
    else /* *p > *q */
        return 1;

}  /* Key_compare */


/********************************************************************/
void Par_bitonic_sort_incr(
        size_t      list_size      /* in     */,
        size_t*    	local_list     /* in/out */,
        size_t*    	local_list1    /* in/out */,
        int       	proc_set_size  /* in     */,
        int 	  	my_rank
        ) {

    unsigned  eor_bit;
    int       proc_set_dim;
    int       stage;
    int       partner;
    //int       my_rank;

    MPI_Comm_rank(COMM_WORLD, &my_rank);

    proc_set_dim = log2(proc_set_size);
    eor_bit = 1 << (proc_set_dim - 1);

    for (stage = 0; stage < proc_set_dim; stage++) {
        partner = my_rank ^ eor_bit;
        if (my_rank < partner){

            Merge_split(
            		list_size,
            		local_list,
            		local_list1,
            		LOW,
            		partner,
            		my_rank
            		);
        }
        else{

            Merge_split(
            		list_size,
            		local_list,
            		local_list1,
            		HIGH,
            		partner,
            		my_rank
            		);
        }
        eor_bit = eor_bit >> 1;
    }
}  /* Par_bitonic_sort_incr */


/********************************************************************/
void Par_bitonic_sort_decr(
        size_t	list_size      /* in     */,
        size_t* local_list     /* in/out */,
        size_t* local_list1    /* in/out */,
        int     proc_set_size  /* in     */,
        int 	my_rank
        ) {

    unsigned  eor_bit;
    int       proc_set_dim;
    int       stage;
    int       partner;
    //int       my_rank;

    MPI_Comm_rank(COMM_WORLD, &my_rank);

    proc_set_dim = log2(proc_set_size);
    eor_bit = 1 << (proc_set_dim - 1);
    for (stage = 0; stage < proc_set_dim; stage++) {
        partner = my_rank ^ eor_bit;

        if (my_rank > partner){
            Merge_split(
            		list_size,
            		local_list,
            		local_list1,
            		LOW,
            		partner,
            		my_rank
            );
        }
        else{
            Merge_split(
            		list_size,
            		local_list,
            		local_list1,
            		HIGH,
            		partner,
            		my_rank
            );
        }
        eor_bit = eor_bit >> 1;
    }

} /* Par_bitonic_sort_decr */


/********************************************************************/
void Merge_split(
        size_t    list_size     /* in     */,
        size_t    *local_list1   /* in/out */,
        size_t	  *local_list2  /* in/out */,
        int       which_keys    /* in     */,
        int       partner       /* in     */,
        int 	  rank			/* in 	  */) {

	 int number_amount;
	 size_t k=0;

	 size_t *temp_key_list1 = malloc(list_size*sizeof(size_t));
	 size_t *temp_key_list2 = malloc(list_size*sizeof(size_t));

	 assert(temp_key_list1 != 0);
	 assert(temp_key_list2 != 0);

	 temp_key_list1 = memset(temp_key_list1, 0, sizeof(size_t)*list_size);
	 temp_key_list2 = memset(temp_key_list2, 0, sizeof(size_t)*list_size);

	 /*
	  * we pack the data into 1 vector called interbuff
	  *
	  */

	 MPI_Status status;
	 int res;
	 size_t *interbuff;
	 res = MPI_Alloc_mem((2*list_size)*sizeof(size_t), MPI_INFO_NULL, &interbuff);
	 assert(res == MPI_SUCCESS);
	 size_t *pos_buff=interbuff;

	 for ( k = 0 ; k < list_size; k++ ){
		 interbuff[k] 				=  local_list1[k];
		 interbuff[k + list_size] 	=  local_list2[k];

	 }

	 size_t *interbuff2;
	 res = MPI_Alloc_mem((2*list_size)*sizeof(size_t), MPI_INFO_NULL, &interbuff2);
	 assert(res == MPI_SUCCESS);


	 MPI_Sendrecv(interbuff,
			 	  2*list_size,
			 	  MPI_LONG_LONG_INT,
			 	  partner,
			 	  0,
			 	  interbuff2,
			 	  2*list_size,
			 	  MPI_LONG_LONG_INT,
			 	  partner,
			 	  0,
			 	  COMM_WORLD,
			 	  &status);


	 MPI_Get_count(&status, MPI_PACKED, &number_amount);

	 for ( k = 0 ; k < list_size; k++ ){
		 temp_key_list1[k] = interbuff2[k];
		 temp_key_list2[k] = interbuff2[k + list_size];

	 }


    if (which_keys == HIGH){
    	Merge_list_high(
    			 list_size,
    			 local_list1,
    			 local_list2,
    			 temp_key_list1,
    			 temp_key_list2
    			 );
    }
    else{
        Merge_list_low(
        		list_size,
        		local_list1,
        		local_list2,
        		temp_key_list1,
        		temp_key_list2
        		);

    }

    if (temp_key_list1) free(temp_key_list1);
    if (temp_key_list2) free(temp_key_list2);

	MPI_Free_mem(interbuff);
	MPI_Free_mem(interbuff2);
    //malloc_trim(0);

} /* Merge_split */


/********************************************************************/
/* Merges the contents of the two lists. */
/* Returns the smaller keys in list1     */
void Merge_list_low(
        size_t   list_size  	/* in     */,
        size_t  *list_key1    	/* in/out */,
        size_t	*list_key2    	/* in/out */,
        size_t  *list_tmp_key1   /* in     */,
        size_t  *list_tmp_key2   /* in     */
        ) {

	size_t  i;
    size_t  index1 = 0;
    size_t  index2 = 0;

    size_t *scratch_list_key1 = malloc(list_size*sizeof(size_t));
    size_t *scratch_list_key2 = malloc(list_size*sizeof(size_t));

    scratch_list_key1[0] = 0;
    scratch_list_key2[0] = 0;

    for (i = 0; i < list_size; i++){
        if (list_key1[index1] <= list_tmp_key1[index2]) {

        	scratch_list_key1[i] = list_key1[index1];
        	scratch_list_key2[i] = list_key2[index1];
        	index1++;

        } else {

        	scratch_list_key1[i] = list_tmp_key1[index2];
        	scratch_list_key2[i] = list_tmp_key2[index2];
        	index2++;
        }
    }
    for (i = 0; i < list_size; i++){
    	list_key1[i] = scratch_list_key1[i];
    	list_key2[i] = scratch_list_key2[i];
    }

    free(scratch_list_key1);
    free(scratch_list_key2);

    //malloc_trim(0);
}  /* Merge_list_low */


/********************************************************************/
/* Returns the larger keys in list 1.    */
void Merge_list_high(
		 size_t   list_size  	/* in     */,
		 size_t  *list_key1    	/* in/out */,
		 size_t	 *list_key2    	/* in/out */,
		 size_t  *list_tmp_key1   /* in     */,
		 size_t	 *list_tmp_key2   /* in     */
		 ) {

    size_t  i;
    size_t  index1 = list_size - 1;
    size_t  index2 = list_size - 1;

    size_t  *scratch_list_key1 = malloc(list_size*sizeof(size_t));
    size_t 	*scratch_list_key2 = malloc(list_size*sizeof(size_t));

    scratch_list_key1[0] =0;
    scratch_list_key2[0]=0;

    size_t counter =0;
    int  rank;
    MPI_Comm_rank(COMM_WORLD, &rank);
    for (i = list_size - 1; i >= 0; i--){

        if (list_key1[index1] >= list_tmp_key1[index2]) {

        	scratch_list_key1[i] = list_key1[index1];
        	scratch_list_key2[i] = list_key2[index1];

        	index1--;
        	counter++;

        } else {

        	scratch_list_key1[i] = list_tmp_key1[index2];
        	scratch_list_key2[i] = list_tmp_key2[index2];

            index2--;
            counter++;
        }

        if (counter >= list_size)
            break;
    }

    for (i = 0; i < list_size; i++){

        	list_key1[i]  = scratch_list_key1[i];
        	list_key2[i]  = scratch_list_key2[i];

    }

    free(scratch_list_key1);
    free(scratch_list_key2);

}  /* Merge_list _high */

/*
 * -------------------------            qksort          ------------------------------
 */

int bitonic_qksort(void *data, size_t size, size_t esize, size_t i, size_t k, int (*compare)(const void *key1, const void *key2)){

	size_t j;

	/*
	 * stop recursion when it is not possible to partition further
	 * when calling qksort:
	 * i = 0
	 * k = size-1
	 */

	while (i < k){

		/*
		 * find where to partition the elements
		 */

		if ((j = bitonic_partition(data, esize, i, k, compare)) < 0){
			return -1;
		}

		/*
		 * recursively sort the left partition
		 */
		if (bitonic_qksort(data, size, esize, i, j, compare) < 0)
			return -1;

		/*
		 * iterate and sort the right partition
		 */
		i = j + 1;
	}
	return 0;
}


int bitonic_issort(void *data, size_t size, size_t esize, int (*compare)(const void *key, const void *key2)){

	size_t *a = data;
	size_t *key;
	size_t i,j;

	if ((key = malloc(sizeof(size_t))) == NULL)
		return -1;

	for ( j = 1; j < size; j++){
		memcpy(key, &a[j], sizeof(size_t));
		i = j - 1;
		while (i >= 0 && compare(&a[i], key) > 0){
			memcpy(&a[(i + 1)], &a[i], sizeof(size_t));
			i--;
		}
		memcpy(&a[(i + 1)], key, sizeof(size_t));
	}

	free(key);
	return 0;

}

int bitonic_partition(void *data, size_t esize, size_t i, size_t k, int (*compare)(const void *key1, const void *key2)){

	size_t *a = data;
	size_t *pval, *temp;
	size_t r[3];
	/*
	 * allocate value for the partition value and swapping
	 */
	if ((pval = malloc(sizeof(size_t))) == NULL)
		return -1;
	if ((temp = malloc(sizeof(size_t))) == NULL){
		free(pval);
		return -1;
	}


	/*
	 * Use the median-of-three method to find partition value
	 */

	int part = (k - i + 1);
	if(part <= 0)
		part = 1;

	r[0] = (rand() % part + i);
	r[1] = (rand() % part + i);
	r[2] = (rand() % part + i);

	/*
	 * TODO: replace the qsort with issort
	 *
	 * issort(r, 3, sizeof(size_t), compare_size_t_V2);
	 */
	qsort(r, 3, sizeof(size_t),compare_size_t_V2);
	memcpy(pval, &a[r[1]], sizeof(size_t));

	/*
	 * Create 2 partitions around the partition value
	 */
	i--;
	k++;
	while(1) {

		/*
		 * move left until an element is found in the wrong partition
		 */

		do {
			k--;
		} while (compare(&a[k], pval) > 0);

		/*
		 * move right until an element is found in the wrong partition
		 */

		do {
			i++;
		} while (compare(&a[i], pval) < 0);

		if (i >= k){
			/*
			 * break when left and right counter cross
			 */
			break;
		}

		else{
			// swap element under the left and right counters
			memcpy(temp, &a[i], sizeof(size_t));
			memcpy(&a[i], &a[k], sizeof(size_t));
			memcpy(&a[k], temp, sizeof(size_t));
		}
	}

	/*
	 * free the storage allocated for partitioning
	 */
	free(pval);
	free(temp);

	/*
	 * return position dividing the two partition
	 */
	return k;

}
