#include <stdlib.h>
#include "tools.c"

static double ** calcWAM(double ** vectors, int n, int vectorDimension) {
	int i, j, d;
	double ** Matrice;
	Matrice = (double**)malloc(n * sizeof(double*));
	if (Matrice == NULL) return NULL;
	d = vectorDimension;
	for (i = 0; i < n; ++i) {
		Matrice[i] = (double*)malloc(n * sizeof(double));
		if (Matrice[i] == NULL) return NULL;
	}
	for (i = 0; i < n; ++i) {
		for (j = i; j < n; ++j) {
			if (j == i) Matrice[i][j] = 0;
			else Matrice[i][j] = Matrice[j][i] =
				exp(-distance(vectors[i], vectors[j], d) / 2.0);
		}
	}
	return Matrice;
}
