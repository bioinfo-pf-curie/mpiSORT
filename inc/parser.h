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
     parser.h

   Authors:
    Frederic Jarlier from Institut Curie
	Nicolas Fedy from Institut Curie
	Leonor Sirotti from Institut Curie
	Thomas Magalhaes from Institut Curie
	Paul Paganiban from Institut Curie
*/


#ifndef PARSER
	#define PARSER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include "tokenizer.h"
#include "time.h"

#define MAX_LINE_SIZE 2048
#define UNMAPPED "unmapped"
#define DISCORDANT "discordant"
typedef struct Flags Flags;
typedef struct Read Read;

#define MODE_NAME	0
#define MODE_OFFSET	1
extern char parse_mode;

struct Flags
{
	unsigned char chr : 5;
	unsigned char is_mate : 1;
	unsigned char left : 1;
	unsigned char replace_gene_with_mgene:1;
};

struct Read
{
	/*
	 * task: REMOVE STRING
	 */
	char* string;
	size_t coord;
	size_t offset_source_file;
	size_t offset;
	char quality;
	struct Read* next;
	struct Read* link;
}; //attribute packed slows by x2 times


typedef struct Read_chain
{
	struct Read* reads;
	struct Read_chain* next;
}Read_chain;

/**
 * \brief Find the header, its size and the chromosomes name
 *
 * \fn unsigned int find_header(char *localData, int rank, size_t *unmappedSize, int *pnbchr, char **pheader, char ***pchrNames)
 * \param localData The local data read directly from the input file
 * \param rank Rank of the process
 * \param unmappedSize Size of lines that are not "Read"
 * \param pnbchr The number of chromosome. Will be set in this function.
 * \param pheader Header in char *
 * \param pchrNames Reference array to the chromosomes names. Will be filled during the function.
 *
 * \return Size of the header, chromosome name
 */
unsigned int find_header(char *localData, int rank, size_t *unmappedSize, int *pnbchr, char **pheader, char ***pchrNames);

/**
 * \brief Initialize the start offset for each process
 *
 * \fn size_t * init_goff(MPI_File mpi_filed,unsigned int headerSize,size_t fsize,int numproc,int rank)
 * \param mpi_filed File to read
 * \param headerSize The size of the header
 * \param fsize The size of the mpi_filed
 * \param numproc The number of process
 * \param rank Rank of the process
 *
 * \return An array of "numproc" cells with the first offset to read for each process (last cell is the last offset of mpi_filed = fsize)
 */
size_t * init_goff(MPI_File mpi_filed,unsigned int headerSize,size_t fsize,int numproc,int rank);

/**
 * \brief Parse localData to reads. This functions takes care of the mate of each read.
 *
 * \fn void parser_paired(char *localData, int rank, size_t start_offset, unsigned char threshold,int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads)
 * \param localData The local data read directly from the input file
 * \param rank Rank of the process
 * \param start_offset First localData offset in the file
 * \param threshold Minimum quality that the reads have to be. If the read quality isn't high enough, it won't be keeped.
 * \param nbchrom The number of chromosome
 * \param preadNumberByChr Reference array to the number of reads by chromosome.
 * \param chrNames Reference array to the chromosomes names.
 * \param preads The reads linked list. Will be set in this function.
 *
 * \return void
 */
void parser_paired(char *localData, int rank, size_t start_offset, unsigned char threshold,int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads);

void parser_single(char *localData, int rank, size_t start_offset, unsigned char threshold,int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads);


/**
 * \brief Extract integer number of current chromosome given as a string
 *
 * \param str The string containing the chromosome ID
 * \param chrNames Reference array to the chromosomes names
 * \param nbchr The number of chromosome
 *
 * \return The chromosome as an integer
 */
int getChr(char* str, char** chrNames, int nbchr);

#endif
