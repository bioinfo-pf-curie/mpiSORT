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
     tokenizer.c

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/



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
