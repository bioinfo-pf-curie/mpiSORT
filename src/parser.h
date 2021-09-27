/*
   This file is part of mpiSORT
   
   Copyright Institut Curie 2020
   
   This software is a computer program whose purpose is to sort SAM file.
   
   You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
   
   The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
   
   The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/

/*
   Module:
     parser.h

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/

#ifndef PARSER
	#define PARSER

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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
 * \brief Initialize the start offset for each process
 *
 * \fn size_t * init_goff(int mpi_file_descriptor,unsigned int headerSize,size_t fsize,int numproc,int rank)
 * \param mpi_file_descriptor File to read
 * \param headerSize The size of the header
 * \param fsize The size of the mpi_file_descriptor
 * \param numproc The number of process
 * \param rank Rank of the process
 *
 * \return An array of "numproc" cells with the first offset to read for each process (last cell is the last offset of mpi_file_descriptor = fsize)
 */
void init_goff(int mpi_file_descriptor,unsigned int headerSize,size_t fsize,int numproc,int rank, size_t *goff);

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
