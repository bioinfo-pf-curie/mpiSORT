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
     merge_utils.c

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/



#include "merge_utils.h"

void indexing(int rank, size_t readNum, Read* chr, size_t *dsend[]){

	size_t j;

	dsend[0] = (size_t*)malloc(readNum*sizeof(size_t));
	dsend[1] = (size_t*)malloc(readNum*sizeof(size_t));
	dsend[2] = (size_t*)malloc(readNum*sizeof(size_t));


	for(j = 0; j < readNum; j++){

		dsend[0][j] = (size_t)rank;
		dsend[1][j] = chr->coord;
		dsend[2][j] = chr->offset;
		chr = chr->next;
	}
}

/**
 * \fn int number_sons(int rank, int num_proc)
 * \brief calculates the number of sons
 * \param rank is the process' rank
 * \param num_proc is the processes total number
 */

int number_sons(int rank, int num_proc){

	int res = 0;

	if (2 * rank + 2 < num_proc){
		res = 2;
	}

	else if (2 * rank + 1 < num_proc){
		res = 1;
	}

	else {
		res = 0;
	}

	return res;
}

void arrayMerge(size_t* a1[], size_t* a2[], size_t* a3[], size_t size1, size_t size2){
	size_t i, j, k;
	i = 0;
	j = 0;
	k = 0;

	while (i < size1 && j < size2){
		if (a1[1][i] < a2[1][j]){
			a3[0][k] = a1[0][i];
			a3[1][k] = a1[1][i];
			a3[2][k] = a1[2][i];


			i++;
			k++;
		}
		else {
			a3[0][k] = a2[0][j];
			a3[1][k] = a2[1][j];
			a3[2][k] = a2[2][j];


			j++;
			k++;
		}
	}

	while (i < size1){
		a3[0][k] = a1[0][i];
		a3[1][k] = a1[1][i];
		a3[2][k] = a1[2][i];


		i++;
		k++;
	}

	while (j < size2){
		a3[0][k] = a2[0][j];
		a3[1][k] = a2[1][j];
		a3[2][k] = a2[2][j];

		j++;
		k++;
	}
}

size_t bucket_arrayMerge(size_t* a1[], size_t* a2[], size_t* a3[], size_t* a4[], size_t s1, size_t s2, size_t s3, char lacking, size_t disp){

	size_t i, j, k, m;
	m = i = j = 0;
	k  = disp;

	while(i < s1 && j < s2 && k < s3){

		if(a1[1][i] < a2[1][j]){
			if(a1[1][i] < a3[1][k]){

				cellcpy(a1, a4, i, m);

				i++;
				m++;
			}

			else{

				cellcpy(a3, a4, k, m);

				k++;
				m++;
			}
		}

		else if(a2[1][j] < a3[1][k]){

			cellcpy(a2, a4, j, m);

			j++;
			m++;
		}

		else{

			cellcpy(a3, a4, k, m);

			k++;
			m++;
		}
	}

	while (i < s1 && j < s2){
		if (a1[1][i] < a2[1][j]){

			cellcpy(a1, a4, i, m);

			i++;
			m++;
		}
		else {

			cellcpy(a2, a4, j, m);

			j++;
			m++;
		}
	}

	while (i < s1 && k < s3){
		if (a1[1][i] < a3[1][k]){

			cellcpy(a1, a4, i, m);

			i++;
			m++;
		}
		else {

			cellcpy(a3, a4, k, m);

			k++;
			m++;
		}
	}

	while (j < s2 && k < s3){
		if (a1[1][j] < a3[1][k]){

			cellcpy(a2, a4, j, m);

			j++;
			m++;
		}
		else {

			cellcpy(a3, a4, k, m);

			k++;
			m++;
		}
	}

	while (i < s1){

		cellcpy(a1, a4, i, m);

		i++;
		m++;
	}

	while (j < s2){

		cellcpy(a2, a4, j, m);

		j++;
		m++;
	}

	if(!lacking){ //if a son has not sent all of its data yet, the following iteration must not occur
		while (k < s3){

			cellcpy(a3, a4, k, m);

			k++;
			m++;
		}
	}

	return k;
}

void cellcpy(size_t* src[], size_t* dst[], size_t i1, size_t i2){

	dst[0][i2] = src[0][i1];
	dst[1][i2] = src[1][i1];
	dst[2][i2] = src[2][i1];
}
