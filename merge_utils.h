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
     merge_utils.h

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
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
