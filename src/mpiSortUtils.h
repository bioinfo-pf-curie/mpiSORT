/*
   This file is part of mpiSORT
   
   Copyright Institut Curie 2022
   
   This software is a computer program whose purpose is to sort SAM file.
   
   You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
   
   The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
   
   The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/

/*
   Module:
     mpiSortUtils.h

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/

#ifndef MPISORTUTILS_H
#define MPISORTUTILS_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <parser.h>

/**************************************************
 *
 *	list interface to handle chromosoms name
 *	collision in the hash table
 *  
 **************************************************/


typedef struct Elmt_ {
		int index; //hold the chromosom index
} Elmt; 

/**************************************************
 *
 *  Hash function prototype for chromosoms string name
 *  
 **************************************************/

#define 	     PRIME_TBLSIZ	65531
unsigned int hashpjw(const void *key);
/**************************************************
 *
 *  Hash function prototype for chromosoms string name
 *  
 **************************************************/


typedef struct CHTbl_ {
	unsigned int (*h)(const void *key);
        int buckets;
        int size;
        Elmt *elements;

} CHTbl;

void print_table(CHTbl *htbl);

int chtbl_init(CHTbl *htbl, int buckets, unsigned int (*h)(const void *key));

void chtbl_destroy(CHTbl *htbl);

int chtbl_insert(CHTbl *htbl, const void *data, int i);

int chtbl_remove(CHTbl *htbl, char *data);

int chtbl_lookup(const CHTbl *htbl, char *data);

#define chtbl_size(htbl) ((htbl)->size)


/**************************************
 *
 *      Other functions of mpiSort.c
 *
 *


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

*/




#endif

