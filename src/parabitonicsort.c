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
     parabitonicsort.c

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/



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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mpi.h"
#include "parabitonicsort.h"

// we limit to 1gb per proc
//#define MAX 1024*1024

#define LOW 0
#define HIGH 1

// the number of key max is 1gb
//#define KEY_MAX 1024*1024
#define key_mpi_t MPI_LONG_LONG_INT
#define index_mpi_t MPI_LONG_LONG_INT

static MPI_Comm COMM_WORLD;

//KEY_T temp_key_list[MAX]; /* buffer for keys received */
                            /* in Merge_split           */

//KEY_T temp_index_list[MAX]; /* buffer for index received */
                      	  	  /* in Merge_split            */

//KEY_T scratch_list_key[MAX]; /* temporary storage for */
                               /* merges                */

//KEY_T scratch_list_index[MAX];


/********************************************************************/

size_t *base_arr2;
size_t rank_to_display = 6;

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

	 if (base_arr2[aa] > base_arr2[bb])
		return 1;
	else if (base_arr2[aa] < base_arr2[bb])
		return -1;
	else
		return 0;
}

void ParallelBitonicSort(MPI_Comm split_comm, int my_rank, int dimension,
		size_t *local_list, size_t *local_index,
			size_t list_size, size_t zero_padding) {

    int       proc_set_size;
    unsigned  and_bit;
    int k = 0;
    COMM_WORLD = split_comm;
    //we add zero at the end of the local_list
    //fprintf(stderr, "Rank %d :::::[BITONIC SORT] list_size  = %zu \n", my_rank, list_size);
    //fprintf(stderr, "Rank %d :::::[BITONIC SORT] zero_padding  = %zu \n", my_rank, zero_padding);

    if ( my_rank == 3 || my_rank == 2 ){
		for (k = 0; k < (zero_padding + 5); k++){
			// fprintf(stderr, "Rank %d :::::[BITONIC SORT] first elements at %zu in local_list = %zu \n", my_rank, k, local_list[k]);
		 }

		for (k = 0; k < (zero_padding + 5); k++){
			// fprintf(stderr, "Rank %d :::::[BITONIC SORT] last elements at %zu in local_list = %zu \n", my_rank, (list_size - k),local_list[list_size - k]);
		}
    }

    if (my_rank < (dimension - 1)){
    	for (k = 0; k < zero_padding; k++){
    		local_list[list_size - k - 1] = 0;
    	}
    }

    /*
    if (my_rank == 1) {
    	for (k = 0; k < 10; k++){
    		fprintf(stderr, "Rank %d :::::[BITONIC SORT] first elements before local_sort at %d local_list = %zu \n", my_rank, k, local_list[k]);
    		fprintf(stderr, "Rank %d :::::[BITONIC SORT] first indexes before local_sort at %d  local_index = %zu \n", my_rank, k, local_index[k]);
    	}
   		for (k = list_size - 1; k > (list_size - 10); k--){
   			fprintf(stderr, "Rank %d :::::[BITONIC SORT] last elements before local_sort at %d local_list = %zu \n", my_rank, k, local_list[k]);
   			fprintf(stderr, "Rank %d :::::[BITONIC SORT] last indexes before local_sort at %d local_index = %zu \n", my_rank, k, local_index[k]);
   		}
    }
	*/
    Local_sort(list_size, local_list, local_index);
    /*
    if (my_rank == 1) {

		for (k = 0; k < 10; k++){
			fprintf(stderr, "Rank %d :::::[BITONIC SORT] first elements after local_sort at %d local_list = %zu \n", my_rank, k, local_list[k]);
			fprintf(stderr, "Rank %d :::::[BITONIC SORT] first indexes after local_sort at %d  local_index = %zu \n", my_rank, k, local_index[k]);
		}


		for (k = list_size - 1; k > (list_size - 10); k--){
			fprintf(stderr, "Rank %d :::::[BITONIC SORT] last elements after local_sort at %d local_list = %zu \n", my_rank, k, local_list[k]);
			fprintf(stderr, "Rank %d :::::[BITONIC SORT] last indexes after local_sort at %d local_index = %zu \n", my_rank, k, local_index[k]);
		}
    }
	*/

    //fprintf(stderr, "Rank %d :::::[BITONIC SORT] after local sort  \n", my_rank);
    /* and_bit is a bitmask that, when "anded" with  */
    /* my_rank, tells us whether we're working on an */
    /* increasing or decreasing list                 */
    for (proc_set_size = 2, and_bit = 2; proc_set_size <= dimension;
    		proc_set_size = proc_set_size*2, and_bit = and_bit << 1){

        if ((my_rank & and_bit) == 0){

            Par_bitonic_sort_incr(list_size, local_list, local_index, proc_set_size);
            //fprintf(stderr, "Rank %d :::::[BITONIC SORT] after Par_bitonic_sort_incr  \n", my_rank);
        }
        else{

            Par_bitonic_sort_decr(list_size, local_list, local_index, proc_set_size);
            //fprintf(stderr, "Rank %d :::::[BITONIC SORT] after Par_bitonic_sort_decr  \n", my_rank);
        }
    }


    //fprintf(stderr, "Rank %d :::::[BITONIC SORT] FINISH BITONIC SORT  \n", my_rank);
    /*
    if (my_rank == 1) {

		for (k = 0; k < 10; k++){
			fprintf(stderr, "Rank %d :::::[BITONIC SORT] first elements after bitonic_sort at %d local_list = %zu \n", my_rank, k, local_list[k]);
			fprintf(stderr, "Rank %d :::::[BITONIC SORT] first indexes after bitonic_sort at %d  local_index = %zu \n", my_rank, k, local_index[k]);
		}

		for (k = list_size - 1; k > (list_size - 10); k--){
			fprintf(stderr, "Rank %d :::::[BITONIC SORT] last elements after bitonic_sort at %d local_list = %zu \n", my_rank, k, local_list[k]);
			fprintf(stderr, "Rank %d :::::[BITONIC SORT] last indexes after bitonic_sort at %d local_index = %zu \n", my_rank, k, local_index[k]);
		}
    }
    */

}


/*********************************************************************/
void Local_sort(
         size_t   list_size     /* in     */,
         size_t  *local_keys  /* in/out */,
         size_t  *local_index /* in/out */) {

	//we create an index vector
	size_t *local_keys_temp = (size_t *)malloc(sizeof(size_t)*list_size);
	size_t *local_index_temp = (size_t *)malloc(sizeof(size_t)*list_size);
	size_t *index_vector = (size_t *)malloc(sizeof(size_t)*list_size);;

	size_t j = 0;

	for(j = 0; j < list_size; j++){
			index_vector[j] = j;
	}

	base_arr2 = local_keys;
	bitonic_qksort(index_vector, list_size, sizeof(size_t), 0, list_size - 1, compare_size_t);

	//then we apply loac index to local_keys
	for(j = 0; j < list_size; j++){
		local_keys_temp[j] = local_keys[index_vector[j]];
		local_index_temp[j] = local_index[index_vector[j]];
	}

	for(j = 0; j < list_size; j++){
		local_keys[j] = local_keys_temp[j];
		local_index[j] = local_index_temp[j];
	}
	free(index_vector);
	free(local_keys_temp);
	free(local_index_temp);
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
int log_base2(int x) {
    int count = 0;

    while (x > 1) {
        x = x/2;
        count++;
    }

    return count;

}  /* log_base2 */


/********************************************************************/
void Par_bitonic_sort_incr(
        int       list_size      /* in     */,
        size_t*    local_list    /* in/out */,
        size_t*	  local_index    /* in/out */,
        int       proc_set_size  /* in     */) {

    unsigned  eor_bit;
    int       proc_set_dim;
    int       stage;
    int       partner;
    int       my_rank;

    MPI_Comm_rank(COMM_WORLD, &my_rank);

    proc_set_dim = log_base2(proc_set_size);
    eor_bit = 1 << (proc_set_dim - 1);

    for (stage = 0; stage < proc_set_dim; stage++) {
        partner = my_rank ^ eor_bit;
        if (my_rank < partner){
            Merge_split(list_size, local_list, local_index, LOW, partner, my_rank);
        }
        else{
            Merge_split(list_size, local_list, local_index, HIGH, partner, my_rank);
        }
        eor_bit = eor_bit >> 1;
    }
}  /* Par_bitonic_sort_incr */


/********************************************************************/
void Par_bitonic_sort_decr(
        int       list_size      /* in     */,
        size_t*    local_list     /* in/out */,
        size_t*	  local_index    /* in/out */,
        int       proc_set_size  /* in     */) {

    unsigned  eor_bit;
    int       proc_set_dim;
    int       stage;
    int       partner;
    int       my_rank;

    MPI_Comm_rank(COMM_WORLD, &my_rank);

    proc_set_dim = log_base2(proc_set_size);
    eor_bit = 1 << (proc_set_dim - 1);
    for (stage = 0; stage < proc_set_dim; stage++) {
        partner = my_rank ^ eor_bit;
        if (my_rank > partner){
            Merge_split(list_size, local_list, local_index, LOW, partner, my_rank);
        }
        else{
            Merge_split(list_size, local_list, local_index, HIGH, partner, my_rank);
        }
        eor_bit = eor_bit >> 1;
    }

} /* Par_bitonic_sort_decr */


/********************************************************************/
void Merge_split(
        size_t       list_size     /* in     */,
        size_t     *local_list  /* in/out */,
        size_t     *local_index /* in/out */,
        int       which_keys    /* in     */,
        int       partner       /* in     */,
        int rank				/* in 	  */) {

	int number_amount;
    MPI_Status status;

    /* key_mpi_t is an MPI (derived) type */
    /* send recieve on local_list*/

    size_t *temp_key_list = (size_t *)malloc(list_size*sizeof(size_t));
    size_t *temp_index_list = (size_t *)malloc(list_size*sizeof(size_t));

    //inititalization
    temp_key_list[0] = 0;
    temp_index_list[0] = 0;

    MPI_Sendrecv(local_list, list_size, MPI_LONG_LONG_INT,
                 partner, 0, temp_key_list, list_size,
                 MPI_LONG_LONG_INT, partner, 0, COMM_WORLD, &status);


    MPI_Sendrecv(local_index, list_size, MPI_LONG_LONG_INT,
                     partner, 0, temp_index_list, list_size,
                     MPI_LONG_LONG_INT, partner, 0, COMM_WORLD, &status);

    MPI_Get_count(&status, MPI_INT, &number_amount);
    //fprintf(stderr, "Rank %d :::::[BITONIC SORT][MERGE SPLIT] After SENDRECV 2 amount = %d \n", partner, number_amount);

    if (which_keys == HIGH){
    	 //fprintf(stderr, "Rank %d :::::[BITONIC SORT][MERGE SPLIT] Before Merge list high \n", rank);
    	 Merge_list_high(list_size, local_list, local_index, temp_key_list, temp_index_list);
    	 //fprintf(stderr, "Rank %d :::::[BITONIC SORT][MERGE SPLIT] After Merge list high \n", rank);
    }
    else{
    	//fprintf(stderr, "Rank %d :::::[BITONIC SORT][MERGE SPLIT] Before Merge list low \n", rank);
        Merge_list_low(list_size, local_list, local_index, temp_key_list, temp_index_list);
        // fprintf(stderr, "Rank %d :::::[BITONIC SORT][MERGE SPLIT] After Merge list low \n", rank);
    }
    free(temp_key_list);
    free(temp_index_list);

} /* Merge_split */


/********************************************************************/
/* Merges the contents of the two lists. */
/* Returns the smaller keys in list1     */
void Merge_list_low(
        size_t   list_size  	/* in     */,
        size_t *list_key    	/* in/out */,
        size_t  *list_index    	/* in/out */,
        size_t  *list_tmp_key   /* in     */,
        size_t  *list_tmp_index /* in     */) {

	size_t  i;
    size_t  index1 = 0;
    size_t  index2 = 0;

    size_t *scratch_list_key = (size_t *)malloc(list_size*sizeof(size_t));
    size_t *scratch_list_index = (size_t *)malloc(list_size*sizeof(size_t));

    for (i = 0; i < list_size; i++){
        if (list_key[index1] <= list_tmp_key[index2]) {

        	scratch_list_key[i] = list_key[index1];
            scratch_list_index[i] = list_index[index1];
            index1++;

        } else {

        	scratch_list_key[i] = list_tmp_key[index2];
            scratch_list_index[i] = list_tmp_index[index2];
            index2++;

        }
    }
    for (i = 0; i < list_size; i++){
    	list_key[i] = scratch_list_key[i];
    	list_index[i] = scratch_list_index[i];
    }

    free(scratch_list_key);
    free(scratch_list_index);
}  /* Merge_list_low */


/********************************************************************/
/* Returns the larger keys in list 1.    */
void Merge_list_high(
		size_t   list_size  /* in     */,
		size_t  *list_key    /* in/out */,
		size_t  *list_index    /* in/out */,
		size_t  *list_tmp_key    /* in     */,
		size_t  *list_tmp_index    /* in     */) {


    size_t  i;
    size_t  index1 = list_size - 1;
    size_t  index2 = list_size - 1;
    size_t *scratch_list_key = (size_t *)malloc(list_size*sizeof(size_t));
    scratch_list_key[0]=0;

    size_t *scratch_list_index = (size_t *)malloc(list_size*sizeof(size_t));
    scratch_list_index[0]=0;

    size_t counter =0;
    int  rank;
    MPI_Comm_rank(COMM_WORLD, &rank);

    /*
    fprintf(stderr, "Rank %d :::::[BITONIC SORT][MERGE LIST HIGH] In Merge list high \n", rank);
    fprintf(stderr, "Rank %d :::::[BITONIC SORT][MERGE LIST HIGH] In Merge list high  list_size = %zu \n", rank, list_size);
    fprintf(stderr, "Rank %d :::::[BITONIC SORT][MERGE LIST HIGH] In Merge list high loop 1 \n", rank);
	*/
    //assert(list_size <= 47756);

    for (i = list_size - 1; i >= 0; i--){

        if (list_key[index1] >= list_tmp_key[index2]) {

        	scratch_list_key[i] = list_key[index1];
        	scratch_list_index[i] = list_index[index1];
        	index1--;

        	//assert(index1 >= 0);
        	/*
        	if (index1 >= list_size){
        	     fprintf(stderr, "Rank %d :::::[BITONIC SORT][MERGE LIST HIGH] PROBLEM In Merge list high index1 = %zu \n", rank, index1);
        	     fprintf(stderr, "Rank %d :::::[BITONIC SORT][MERGE LIST HIGH] PROBLEM In Merge list high counter = %zu \n", rank, counter);
        	}
        	*/
        	//assert(index1 < list_size);
        	counter++;

        } else {

        	scratch_list_key[i] = list_tmp_key[index2];
        	scratch_list_index[i] = list_tmp_index[index2];
            index2--;
            /*
            if (index2 >= list_size){
                  fprintf(stderr, "Rank %d :::::[BITONIC SORT][MERGE LIST HIGH] PROBLEM In Merge list high index2 = %zu \n", rank, index1);
                  fprintf(stderr, "Rank %d :::::[BITONIC SORT][MERGE LIST HIGH] PROBLEM In Merge list high counter = %zu \n", rank, counter);
             }
			*/
            counter++;
        }

        if (counter >= list_size)
            break;
    }

    //fprintf(stderr, "Rank %d :::::[BITONIC SORT][MERGE LIST HIGH] In Merge list high loop 2 \n", rank);
    for (i = 0; i < list_size; i++){

        	list_key[i] = scratch_list_key[i];
        	list_index[i] = scratch_list_index[i];

    }

    free(scratch_list_key);
    free(scratch_list_index);
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

	if ((key = (size_t *)malloc(sizeof(size_t))) == NULL)
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
	if ((pval = (size_t *)malloc(sizeof(size_t))) == NULL)
		return -1;
	if ((temp = (size_t *)malloc(sizeof(size_t))) == NULL){
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
