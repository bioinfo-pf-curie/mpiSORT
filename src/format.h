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

#ifndef FORMAT_H_
	#define FORMAT_H_

	size_t** format(int rank, int num_proc, size_t* recv[3], size_t readNum, size_t** count_diffusep, size_t *abs_offset);
	void countNumRead(int num_proc, size_t* recv[3], size_t readNum, size_t* count_diffuse);

#endif
