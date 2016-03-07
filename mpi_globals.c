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
     mpi_globals.c

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/



#include <errno.h>
#include <stdio.h>
#include <limits.h>
#include "mpi_globals.h"
#include "mpi.h"

int g_rank;
int g_size;
int g_iter = 0;

char g_name[MPI_MAX_PROCESSOR_NAME];

FILE * g_fd = NULL;

double g_start;
double g_end;
double g_elapsed = 0; // time elapised measurement

void mpi_print_free_memory(void)
{
	FILE * fp = popen("free -m | grep Mem | awk '{print $2}'", "r");
	if (NULL == fp)
	{
		// errno 12 means not enough memory
		DEBUG("COULD NOT READ FREE MEMORY, error=%d\n", errno);
		return;
	}
	char path[1035];
	while (fgets(path, sizeof(path)-1, fp) != NULL) {
		DEBUG(" \n Free MB: %s \n", path);
	}

	DEBUG("\n its done!! \n");
	pclose(fp);
}



int MPI_Large_bcast(void * buffer, unsigned long count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
	while (count > INT_MAX)
	{
		if (MPI_SUCCESS != MPI_Bcast(buffer, INT_MAX, datatype, root, comm))
		{
			DEBUG("Failed to broadcast %d count\n", INT_MAX);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		count -= INT_MAX;
		buffer += INT_MAX;
	}

	if (MPI_SUCCESS != MPI_Bcast(buffer, count, datatype, root, comm))
	{
		DEBUG("Failed to broadcast %lu count\n", count);
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	return MPI_SUCCESS;
}

void mpi_open_psam_file()
{
	char filename[32];

	snprintf(&filename[0], 32, "psam/iter%dproc%d.pSAM", g_iter, g_rank);

	if (NULL != g_fd)
	{
		DEBUG("Global file descriptor was not NULL; cannot continue.\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	if (NULL == (g_fd = fopen(filename, "w")))
	{
		DEBUG("Failed to open file: %s\n", filename);
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
}

void mpi_close_psam_file()
{
	if (0 != fclose(g_fd))
	{
		DEBUG("Failed to close file descriptor\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	g_fd = NULL;
}
