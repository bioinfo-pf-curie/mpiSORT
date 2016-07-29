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
     mpi_globals.h

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/




#ifndef _MPI_GLOBALS_H_
#define _MPI_GLOBALS_H_

#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <errno.h>
#include <mpi.h>
#include <assert.h>

extern int g_rank; // individual processor rank in MPI_COMM_WORLD
extern int g_size; // the total number of processors in MPI_COMM_WORLD
extern int g_iter; // what iterator are we on (used for writing files by each process)
extern FILE * g_fd; // global file descriptor used to write partial SAM file by each process
extern char g_name[MPI_MAX_PROCESSOR_NAME]; // name of node used, ie 'node9'
extern double g_start;
extern double g_end;
extern double g_elapsed;

#define DEBUG(format, ...) fprintf(stderr, "DEBUG %s:%d" format "\n", __FILE__, __LINE__, ##__VA_ARGS__);
#endif

#define FREE_IF_NOT_NULL(x) if((x) != NULL) free (x)

/*
 	char tempString[128];
	char timestamp[100];
	struct timeval tv;
	struct tm *timeptr;
	gettimeofday(&tv, NULL);
	timeptr = localtime(&tv.tv_sec);

	strftime(timestamp, sizeof(timestamp), "%H:%M:%S", timeptr);

	int tempLength = snprintf(&tempString[0], 64, "[%s:%03d - P%d (%s) - %s] ", timestamp,

			(int)(tv.tv_usec / 1000), g_rank, g_name, __FUNCTION__);

	tempLength += snprintf(&tempString[tempLength], 64, format,__VA_ARGS__);
*/

extern void mpi_print_free_memory();
extern void mpi_open_psam_file();
extern void mpi_close_psam_file();
