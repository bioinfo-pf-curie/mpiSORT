/*
   mpiSORT
   Copyright (C) 2016-2017 Institut Curie / Institut Pasteur

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
     tokenizer.c

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

#include <stdio.h>

#include "tokenizer.h"

// does NOT allocate memory : had to be done elsewhere
int tokenizer(char* str, const char delim, char* token){

	static char* last;
	int found = 0, length = 0;
	char* check;
	int i;

	if(str)
		last = str;

	else if (!last || *last == 0)
		return 0;

	check = last;

	while (check && !found && last[length]){

		if (*check == delim)
			found = 1;
		else{
			check++;
			length++;
		}
	}

	if (!found)
		return 0;

	for(i = 0; i < length; i++){
		token[i] = *last++;
	}

	token[length] = 0;
	last++;

	return 1;
}
