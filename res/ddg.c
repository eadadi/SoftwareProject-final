#include <math.h>
#include <stdlib.h>
double **calcDDG(double ** M, int n) {
	double ** Matrice = (double**)calloc(n, sizeof(double*));
	for (int i = 0; i < n; ++i) {
		Matrice[i] = (double*)calloc(n, sizeof(double));
	}

	for (int i = 0; i < n; ++i) {
		double dii = 0;
		for (int j = 0; j < n; ++j) {
			dii += M[i][j];
		}
		Matrice[i][i] = dii;
	}

	return Matrice;
}