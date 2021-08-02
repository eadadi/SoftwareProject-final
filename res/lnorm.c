#include <stdlib.h>
static double ** raise_D_diagonal_in_minus_half(double **D, int n) {
	int i;
	double **raised;
	raised = (double**)calloc(n, sizeof(double*));
	if (raised == NULL) return NULL;
	for (i = 0; i < n; ++i) {
		raised[i] = (double*)calloc(n, sizeof(double));
		if (raised[i] == NULL) return NULL;
	}
	for (i = 0; i < n; ++i) raised[i][i] = pow(D[i][i], -0.5);
	return raised;
}
static double ** left_multip_of_diagonal_matrice(double ** D, double **M, int n) {
	int i, j;
	double **result;
	result = (double**)calloc(n, sizeof(double*));
	if (result == NULL) return NULL;
	for (i = 0; i < n; ++i) {
		result[i] = (double*)calloc(n, sizeof(double));
		if (result[i] == NULL) return NULL;
	}
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			result[i][j] = D[i][i] * M[i][j];
		}
	}
	return result;
}
static double ** right_multip_of_diagonal_matrice(double ** M, double **D, int n) {
	int i, j;
	double **result;
	result = (double**)calloc(n, sizeof(double*));
	if (result == NULL) return NULL;
	for (i = 0; i < n; ++i) {
		result[i] = (double*)calloc(n, sizeof(double));
		if (result[i] == NULL) return NULL;
	}
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			result[i][j] = M[i][j] * D[j][j];
		}
	}
	return result;
}
/* Normalized Graph Laplacian */
static double ** calcNGL(double ** W, double ** D, int n) { 
	int i, j;
	double ** result, **Dtag, **tmp;
	Dtag = raise_D_diagonal_in_minus_half(D, n);
	if (Dtag == NULL) return NULL;
	tmp = left_multip_of_diagonal_matrice(Dtag, W, n);
	if (tmp == NULL) return NULL;
	result = right_multip_of_diagonal_matrice(tmp, Dtag, n);
	if (result == NULL) return NULL;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			if (i == j) result[i][i] = 1 - result[i][i];
			else result[i][j] = -result[i][j];
		}
	}

	for (i = 0; i < n; ++i) free(Dtag[i]);
	free(Dtag);
	for (i = 0; i < n; ++i) free(tmp[i]);
	free(tmp);

	return result;
}