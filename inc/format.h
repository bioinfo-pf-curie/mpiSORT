/*
This file is part of mpiSORT

For sorting NGS data using MPI

The project was developped by

Frederic Jarlier from Institut Curie
Nicolas Fedy from Institut Curie
Leonor Sirotti from Institut Curie
Thomas Magalahes from Institut Curie

Copyright (C) 2016  Institut Curie

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FORMAT_H_
	#define FORMAT_H_

	size_t** format(int rank, int num_proc, size_t* recv[3], size_t readNum, size_t** count_diffusep, size_t *abs_offset);
	void countNumRead(int num_proc, size_t* recv[3], size_t readNum, size_t* count_diffuse);

#endif
