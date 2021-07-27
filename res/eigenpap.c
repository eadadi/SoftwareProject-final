#include <stdlib.h>
#include "value_vector_map_struct.h"


static int cmpfunc(const void *a, const void *b) {

	if (*((value_vector_map*)a)->eigenvalue > *((value_vector_map*)b)->eigenvalue)
		return 1;
	else return -1;
}

static void sortMap(value_vector_map *map, int n) {
	qsort(map, n, sizeof(value_vector_map), cmpfunc);
}

static int determineK(value_vector_map *map, int n) {
	int i, N, k;
	double tmp, max;

	N = n / 2;
	k = 0;
	max = 0;
	for (i = 0; i <= N; ++i) {
		tmp = *map[i + 1].eigenvalue - *map[i].eigenvalue;
		if (tmp > max) {
			tmp = max;
			k = i;
		}
	}
	return k;
}
static value_vector_map* setMap(double ** eigenvectors, double * eigenvalues, int n) {
	int i;
	value_vector_map * map;
	map = (value_vector_map*)calloc(n, sizeof(value_vector_map));
	if (map == NULL) return NULL;
	for (i = 0; i < n; ++i) {
		map[i].eigenvector = &eigenvectors[i];
		map[i].eigenvalue = &eigenvalues[i];
	}
	sortMap(map, n);
	return map;
}
