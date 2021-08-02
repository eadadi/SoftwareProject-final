#include <stdlib.h>
static double **calcDDG(double ** M, int n) {
	double ** Matrice, dii;
	int i, j;
	Matrice = (double**)calloc(n, sizeof(double*));
	if (Matrice == NULL) return NULL;
	for (i = 0; i < n; ++i) {
		Matrice[i] = (double*)calloc(n, sizeof(double));
		if (Matrice[i] == NULL) return NULL;
	}

	for (i = 0; i < n; ++i) {
		dii = 0;
		for (j = 0; j < n; ++j) {
			dii += M[i][j];
		}
		Matrice[i][i] = dii;
	}
	return Matrice;
}
