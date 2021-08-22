#include "../spkmeans.h"

int cmpfunc(const void *a, const void *b) {

	if (*((value_vector_map*)a)->eigenvalue > *((value_vector_map*)b)->eigenvalue)
		return 1;
	else if (*((value_vector_map*)a)->eigenvalue < *((value_vector_map*)b)->eigenvalue)
		return -1;
	else if ((*(value_vector_map*)a).index >= (*(value_vector_map*)b).index)
		return 1;
	else return -1;
}

void sortMap(value_vector_map *map, int n) {
	qsort(map, n, sizeof(value_vector_map), cmpfunc);
}

int determineK(value_vector_map *map, int n) {
	int i, N, k;
	double max;
	double lambda_i, lambda_ipp;
	double delta_i;

	N = n / 2;
	max = 0;
	k = i = 0;
	for (; i <= N; ++i) {
		lambda_i = *(map[i].eigenvalue);
		lambda_ipp = *(map[i + 1].eigenvalue);
		delta_i = fabs(lambda_i - lambda_ipp);
		if (delta_i > max) {
			k = i;
			max = delta_i;
		}
	}

	return k+1;
}


value_vector_map* setMap(double ** eigenvectors, double * eigenvalues, int n) {
	int i;
	value_vector_map * map;
	map = (value_vector_map*)calloc(n, sizeof(value_vector_map));
	if (map == NULL) return NULL;
	for (i = 0; i < n; ++i) {
		map[i].eigenvector = &eigenvectors[i];
		map[i].eigenvalue = &eigenvalues[i];
		map[i].index = i;
	}

	sortMap(map, n);
	return map;
}
