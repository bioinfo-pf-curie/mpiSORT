/*
   mpiSORT
   Copyright (C) 2016-2017 Institut Curie, 26 rue d'Ulm, Paris, France

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
     bufferized_read.c

   Authors:
 	Frederic Jarlier from Institut Curie
	Nicolas Joly from Institut Pasteur
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/

#ifndef QKSORT_H
#define QKSORT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern size_t *base_arr2;

int compare_size_t(const void *a, const void *b);

int compare_size_t_V2(const void *a, const void *b);

int partition(void *data, size_t esize, size_t i, size_t k, int (*compare)(const void *key1, const void *key2));

int qksort(void *data, size_t size, size_t esize, size_t i, size_t k, int (*compare)(const void *key1, const void *key2));

#endif
