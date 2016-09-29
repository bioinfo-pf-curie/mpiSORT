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
