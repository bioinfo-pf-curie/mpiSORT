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
     mergeUtils.h

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/



#ifndef MERGE_UTILS_H_
#define MERGE_UTILS_H_

#include <stdio.h>
#include <stdlib.h>

#include "parser.h"
#include "merge.h"

int number_sons(int rank, int num_proc);
void arrayMerge(size_t* a1[], size_t* a2[], size_t* a3[], size_t size1, size_t size2);
size_t bucket_arrayMerge(size_t* a1[], size_t* a2[], size_t* a3[], size_t* a4[], size_t s1, size_t s2, size_t s3, char lacking, size_t disp);
void cellcpy(size_t* src[], size_t* dst[], size_t i1, size_t i2);

void indexing(int rank, size_t readNum, Read* chr, size_t *send[]);

#endif
