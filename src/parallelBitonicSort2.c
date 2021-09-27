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
     parallelBitonicSort.c

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
#include "malloc.h"
#include "parallelBitonicSort2.h"
#include "mergeSort.h"
// we limit to 1gb per proc

#define LOW 0
#define HIGH 1

// the number of key max is 1gb
#define key_mpi_t MPI_LONG_LONG_INT
#define index_mpi_t MPI_LONG_LONG_INT

static MPI_Comm COMM_WORLD;

/********************************************************************/



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


void ParallelBitonicSort2(
		MPI_Comm split_comm,
		int my_rank,
		int dimension,
		size_t *local_list, 	//stand for coordinates
		int *local_list1, 		//stand for read sizes
		int *local_list2, 		//stand for read rank
		size_t *local_list3,	//stand for read offset source
		int *local_list4,    	//stand for original rank offset source
		size_t list_size) {

	// In this version of bitonic there is no more index.
	// Everything is sorted according
	// the local_list vector.

    int       proc_set_size;
    unsigned  and_bit;
    size_t k = 0;
    COMM_WORLD = split_comm;

    	
    Local_sort2(list_size, local_list, local_list1, local_list2, local_list3, local_list4);

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

            Par_bitonic_sort_incr2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
            		proc_set_size,
            		my_rank);
        }
        else{

            Par_bitonic_sort_decr2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
            		proc_set_size,
            		my_rank);
        }
    }

    for (k = 0; k < (list_size - 1); k++){
      	assert(local_list[k] <= local_list[k + 1]);
    }
}


/*********************************************************************/
void Local_sort2(
         size_t   list_size     /* in     */,
         size_t  *local_keys    /* in/out */,
         int     *local_keys1   /* in/out */,
         int     *local_keys2   /* in/out */,
         size_t  *local_keys3   /* in/out */,
         int     *local_keys4   /* in/out */
         ) {

	//we create an index vector
	size_t *local_keys_temp   	 = malloc(sizeof(size_t)*list_size);
	int    *local_keys_temp1     = malloc(sizeof(int)*list_size);
	int    *local_keys_temp2     = malloc(sizeof(int)*list_size);
	size_t *local_keys_temp3  	 = malloc(sizeof(size_t)*list_size);
	int    *local_keys_temp4  	 = malloc(sizeof(int)*list_size);

	assert(local_keys_temp);
	assert(local_keys_temp1);
	assert(local_keys_temp2);
	assert(local_keys_temp3);
	assert(local_keys_temp4);

	size_t *index_vector = (size_t *)malloc(sizeof(size_t)*list_size);

	size_t j = 0;

	for(j = 0; j < list_size; j++){
			index_vector[j] = j;
	}

	base_arr2 = local_keys;
	//bitonic_qksort2(index_vector, list_size, sizeof(size_t), 0, list_size - 1, compare_size_t);
	MergeSortMain(index_vector, list_size);
	//then we apply loac index to local_keys
	for(j = 0; j < list_size; j++){
		local_keys_temp[j]  = local_keys[index_vector[j]];
		local_keys_temp1[j] = local_keys1[index_vector[j]];
		local_keys_temp2[j] = local_keys2[index_vector[j]];
		local_keys_temp3[j] = local_keys3[index_vector[j]];
		local_keys_temp4[j] = local_keys4[index_vector[j]];
	}

	for(j = 0; j < list_size; j++){
		local_keys[j]  = local_keys_temp[j];
		local_keys1[j] = local_keys_temp1[j];
		local_keys2[j] = local_keys_temp2[j];
		local_keys3[j] = local_keys_temp3[j];
		local_keys4[j] = local_keys_temp4[j];
	}

	free(index_vector);
	free(local_keys_temp);
	free(local_keys_temp1);
	free(local_keys_temp2);
	free(local_keys_temp3);
	free(local_keys_temp4);
	malloc_trim(0);
}


/*********************************************************************/
int Key_compare2(const size_t* p, const size_t* q) {

    if (*p < *q)
        return -1;
    else if (*p == *q)
        return 0;
    else /* *p > *q */
        return 1;

}  /* Key_compare */


/********************************************************************/
void Par_bitonic_sort_incr2(
        size_t      list_size      /* in     */,
        size_t*    	local_list     /* in/out */,
        int*    	local_list1    /* in/out */,
        int*    	local_list2    /* in/out */,
        size_t*    	local_list3    /* in/out */,
        int*    	local_list4    /* in/out */,
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

            Merge_split2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
            		LOW,
            		partner,
            		my_rank
            		);
        }
        else{

            Merge_split2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
            		HIGH,
            		partner,
            		my_rank
            		);
        }
        eor_bit = eor_bit >> 1;
    }
}  /* Par_bitonic_sort_incr */


/********************************************************************/
void Par_bitonic_sort_decr2(
        size_t	list_size      /* in     */,
        size_t* local_list     /* in/out */,
        int*    local_list1    /* in/out */,
        int*    local_list2    /* in/out */,
        size_t* local_list3    /* in/out */,
        int*    local_list4    /* in/out */,
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
            Merge_split2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
            		LOW,
            		partner,
            		my_rank
            );
        }
        else{
            Merge_split2(
            		list_size,
            		local_list,
            		local_list1,
            		local_list2,
            		local_list3,
            		local_list4,
            		HIGH,
            		partner,
            		my_rank
            );
        }
        eor_bit = eor_bit >> 1;
    }

} /* Par_bitonic_sort_decr */


/********************************************************************/
void Merge_split2(
        size_t    list_size     /* in     */,
        size_t    *local_list1   /* in/out */,
        int	      *local_list2  /* in/out */,
        int       *local_list3  /* in/out */,
        size_t    *local_list4  /* in/out */,
        int       *local_list5  /* in/out */,
        int       which_keys    /* in     */,
        int       partner       /* in     */,
        int 	  rank			/* in 	  */) {

	 int number_amount;
	 size_t k=0;

	 size_t *temp_key_list1 	= malloc(list_size*sizeof(size_t));
	 size_t *temp_key_list4   	= malloc(list_size*sizeof(size_t));
	 int *temp_key_list2   		= malloc(list_size*sizeof(int));
	 int *temp_key_list3 		= malloc(list_size*sizeof(int));
	 int *temp_key_list5   		= malloc(list_size*sizeof(int));

	 assert(temp_key_list1 != 0);
	 assert(temp_key_list2 != 0);
	 assert(temp_key_list3 != 0);
	 assert(temp_key_list4 != 0);
	 assert(temp_key_list5 != 0);

	 temp_key_list1 = memset(temp_key_list1, 0, sizeof(size_t)*list_size);
	 temp_key_list2 = memset(temp_key_list2, 0, sizeof(int)*list_size);
	 temp_key_list3 = memset(temp_key_list3, 0, sizeof(int)*list_size);
	 temp_key_list4 = memset(temp_key_list4, 0, sizeof(size_t)*list_size);
	 temp_key_list5 = memset(temp_key_list5, 0, sizeof(int)*list_size);


	 /*
	  * we pack the data into 1 vector called interbuff
	  *
	  */

	 MPI_Status status;
	 int res;
	 size_t *interbuff;
	 res = MPI_Alloc_mem((5*list_size)*sizeof(size_t), MPI_INFO_NULL, &interbuff);
	 assert(res == MPI_SUCCESS);
	 size_t *pos_buff=interbuff;

	 for ( k = 0 ; k < list_size; k++ ){
		 interbuff[k] 				=  local_list1[k];
		 interbuff[k + list_size] 	=  (size_t)local_list2[k];
		 interbuff[k + 2*list_size] =  (size_t)local_list3[k];
		 interbuff[k + 3*list_size] =  local_list4[k];
		 interbuff[k + 4*list_size] =  (size_t)local_list5[k];
	 }

	 size_t *interbuff2;
	 res = MPI_Alloc_mem((5*list_size)*sizeof(size_t), MPI_INFO_NULL, &interbuff2);
	 assert(res == MPI_SUCCESS);


	 MPI_Sendrecv(interbuff,
			 	  5*list_size,
			 	  MPI_LONG_LONG_INT,
			 	  partner,
			 	  0,
			 	  interbuff2,
			 	  5*list_size,
			 	  MPI_LONG_LONG_INT,
			 	  partner,
			 	  0,
			 	  COMM_WORLD,
			 	  &status);


	 MPI_Get_count(&status, MPI_PACKED, &number_amount);

	 for ( k = 0 ; k < list_size; k++ ){
		 temp_key_list1[k] = (size_t)interbuff2[k];
		 temp_key_list2[k] = (int)   interbuff2[k + list_size];
		 temp_key_list3[k] = (int)   interbuff2[k + 2*list_size];
		 temp_key_list4[k] = (size_t)interbuff2[k + 3*list_size];
		 temp_key_list5[k] = (int)   interbuff2[k + 4*list_size];
	 }


    if (which_keys == HIGH){
    	Merge_list_high2(
    			 list_size,
    			 local_list1,
    			 local_list2,
    			 local_list3,
    			 local_list4,
    			 local_list5,
    			 temp_key_list1,
    			 temp_key_list2,
    			 temp_key_list3,
    			 temp_key_list4,
    			 temp_key_list5
    			 );
    }
    else{
        Merge_list_low2(
        		list_size,
        		local_list1,
        		local_list2,
        		local_list3,
        		local_list4,
        		local_list5,
        		temp_key_list1,
        		temp_key_list2,
        		temp_key_list3,
        		temp_key_list4,
        		temp_key_list5
        		);

    }

    if (temp_key_list1) free(temp_key_list1);
    if (temp_key_list2) free(temp_key_list2);
    if (temp_key_list3) free(temp_key_list3);
    if (temp_key_list4) free(temp_key_list4);
    if (temp_key_list5) free(temp_key_list5);

	MPI_Free_mem(interbuff);
	MPI_Free_mem(interbuff2);
    //malloc_trim(0);

} /* Merge_split */


/********************************************************************/
/* Merges the contents of the two lists. */
/* Returns the smaller keys in list1     */
void Merge_list_low2(
        size_t   list_size  	/* in     */,
        size_t  *list_key    	/* in/out */,
        int		*list_key1    	/* in/out */,
        int		*list_key2    	/* in/out */,
        size_t  *list_key3    	/* in/out */,
        int     *list_key4    	/* in/out */,
        size_t  *list_tmp_key   /* in     */,
        int		*list_tmp_key1   /* in     */,
        int		*list_tmp_key2   /* in     */,
        size_t  *list_tmp_key3   /* in     */,
        int     *list_tmp_key4   /* in     */
        ) {

	size_t  i;
    size_t  index1 = 0;
    size_t  index2 = 0;

    size_t *scratch_list_key  = malloc(list_size*sizeof(size_t));
    int	   *scratch_list_key1 = malloc(list_size*sizeof(size_t));
    int    *scratch_list_key2 = malloc(list_size*sizeof(size_t));
    size_t *scratch_list_key3 = malloc(list_size*sizeof(size_t));
    int    *scratch_list_key4 = malloc(list_size*sizeof(int));

    scratch_list_key[0]  = 0;
    scratch_list_key1[0] = 0;
    scratch_list_key2[0] = 0;
    scratch_list_key3[0] = 0;
    scratch_list_key4[0] = 0;

    for (i = 0; i < list_size; i++){
        if (list_key[index1] <= list_tmp_key[index2]) {

        	scratch_list_key[i]  = list_key[index1];
        	scratch_list_key1[i] = list_key1[index1];
        	scratch_list_key2[i] = list_key2[index1];
        	scratch_list_key3[i] = list_key3[index1];
        	scratch_list_key4[i] = list_key4[index1];
            index1++;

        } else {

        	scratch_list_key[i]  = list_tmp_key[index2];
        	scratch_list_key1[i] = list_tmp_key1[index2];
        	scratch_list_key2[i] = list_tmp_key2[index2];
        	scratch_list_key3[i] = list_tmp_key3[index2];
        	scratch_list_key4[i] = list_tmp_key4[index2];
            index2++;
        }
    }
    for (i = 0; i < list_size; i++){
    	list_key[i] = scratch_list_key[i];
    	list_key1[i] = scratch_list_key1[i];
    	list_key2[i] = scratch_list_key2[i];
    	list_key3[i] = scratch_list_key3[i];
    	list_key4[i] = scratch_list_key4[i];

    }

    free(scratch_list_key);
    free(scratch_list_key1);
    free(scratch_list_key2);
    free(scratch_list_key3);
    free(scratch_list_key4);

    //malloc_trim(0);
}  /* Merge_list_low */


/********************************************************************/
/* Returns the larger keys in list 1.    */
void Merge_list_high2(
		 size_t   list_size  	/* in     */,
		 size_t  *list_key    	/* in/out */,
		 int	 *list_key1    	/* in/out */,
		 int	 *list_key2    	/* in/out */,
		 size_t  *list_key3    	/* in/out */,
		 int     *list_key4    	/* in/out */,
		 size_t  *list_tmp_key   /* in     */,
		 int	 *list_tmp_key1   /* in     */,
		 int     *list_tmp_key2   /* in     */,
		 size_t  *list_tmp_key3   /* in     */,
		 int     *list_tmp_key4   /* in     */
		 ) {

    size_t  i;
    size_t  index1 = list_size - 1;
    size_t  index2 = list_size - 1;

    size_t  *scratch_list_key  = malloc(list_size*sizeof(size_t));
    int 	*scratch_list_key1 = malloc(list_size*sizeof(int));
    int 	*scratch_list_key2 = malloc(list_size*sizeof(int));
    size_t  *scratch_list_key3 = malloc(list_size*sizeof(size_t));
    int 	*scratch_list_key4 = malloc(list_size*sizeof(int));

    scratch_list_key[0] =0;
    scratch_list_key1[0]=0;
    scratch_list_key2[0]=0;
    scratch_list_key3[0]=0;
    scratch_list_key4[0]=0;

    size_t counter =0;
    int  rank;
    MPI_Comm_rank(COMM_WORLD, &rank);
    for (i = list_size - 1; i >= 0; i--){

        if (list_key[index1] >= list_tmp_key[index2]) {

        	scratch_list_key[i]  = list_key[index1];
        	scratch_list_key1[i] = list_key1[index1];
        	scratch_list_key2[i] = list_key2[index1];
        	scratch_list_key3[i] = list_key3[index1];
        	scratch_list_key4[i] = list_key4[index1];
        	index1--;
        	counter++;

        } else {

        	scratch_list_key[i]  = list_tmp_key[index2];
        	scratch_list_key1[i] = list_tmp_key1[index2];
        	scratch_list_key2[i] = list_tmp_key2[index2];
        	scratch_list_key3[i] = list_tmp_key3[index2];
        	scratch_list_key4[i] = list_tmp_key4[index2];
            index2--;
            counter++;
        }

        if (counter >= list_size)
            break;
    }

    for (i = 0; i < list_size; i++){

        	list_key[i]   = scratch_list_key[i];
        	list_key1[i]  = scratch_list_key1[i];
        	list_key2[i]  = scratch_list_key2[i];
        	list_key3[i]  = scratch_list_key3[i];
        	list_key4[i]  = scratch_list_key4[i];
    }

    free(scratch_list_key);
    free(scratch_list_key1);
    free(scratch_list_key2);
    free(scratch_list_key3);
    free(scratch_list_key4);
}  /* Merge_list _high */

/*
 * -------------------------            qksort          ------------------------------
 */

int bitonic_qksort2(void *data, size_t size, size_t esize, size_t i, size_t k, int (*compare)(const void *key1, const void *key2)){

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

		if ((j = bitonic_partition2(data, esize, i, k, compare)) < 0){
			return -1;
		}

		/*
		 * recursively sort the left partition
		 */
		if (bitonic_qksort2(data, size, esize, i, j, compare) < 0)
			return -1;

		/*
		 * iterate and sort the right partition
		 */
		i = j + 1;
	}
	return 0;
}


int bitonic_issort2(void *data, size_t size, size_t esize, int (*compare)(const void *key, const void *key2)){

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

int bitonic_partition2(void *data, size_t esize, size_t i, size_t k, int (*compare)(const void *key1, const void *key2)){

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
