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
     mergeSort.c

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/



#include "mergeSort.h"

Read* mergeSort(Read* c, size_t n){
	size_t q, p;
	Read* d;

	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSort(c, p);
		if(q > 1)
			mergeSort(d, q);
	}
	else
		d = c->next;
	d = structMerge(c, p, d, q);

	return d;
}

Read* structMerge(Read* c, size_t p, Read* d, size_t q){

	Read* t;

	while(1){

		if (c->next->coord > d->next->coord){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}
			--p;
		}
		c = c->next;
	}
	return d;
}
