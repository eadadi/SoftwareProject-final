#include <stdlib.h>
#include "tools.c"
static double ** calcWAM(double ** vectors, int n, int vectorDimension) {
	int d = vectorDimension;
	double ** Matrice = (double**)calloc(n, sizeof(double*));
	if (Matrice == NULL) return NULL;
	for (int i = 0; i < n; ++i) {
		Matrice[i] = (double*)calloc(n, sizeof(double));
		if (Matrice[i] == NULL) return NULL;
	}
	for (int i = 0; i < n; ++i) {
		for (int j = i; j < n; ++j) {
			if (j == i) Matrice[i][j] = 0;
			else Matrice[i][j] = Matrice[j][i] =
				exp(-distance(vectors[i], vectors[j], d) / 2);
		}
	}
	return Matrice;
}
