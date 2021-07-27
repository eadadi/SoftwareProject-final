#include<stdlib.h>
#include "value_vector_map_struct.h"
#include "tools.c"
static double ** calcNormalaizedRnk(value_vector_map *map, int n, int clustersNumber) {
	double **Rnk, dist, *zero;
	int i, j;
	Rnk = (double**)malloc(n * sizeof(double*));
	if (Rnk == NULL) return NULL;
	//Obtain first "clusterNumber" eigenvectors
	for (i = 0; i < n; ++i) {
		Rnk[i] = (double*)malloc(clustersNumber * sizeof(double));
		if (Rnk[i] == NULL) return NULL;
	}

	for (j = 0; j < clustersNumber; ++j) {
		for (i = 0; i < n; ++i) {
			Rnk[i][j] = *(*(map[j].eigenvector) + i);
		}
	}
	//Normalize each row
	zero = (double*)calloc(clustersNumber, sizeof(double));
	for (i = 0; i < n; ++i) {
		dist = distance(Rnk[i], zero, clustersNumber);
		for (j = 0; j < clustersNumber; ++j) {
			Rnk[i][j] != 0 ? Rnk[i][j] /= dist : free(NULL);
		}
	}

	free(zero);
	return Rnk;
}