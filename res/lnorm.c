#include "../spkmeans.h"
void raise_D_diagonal_in_minus_half(double **D, int n) {
	int i;
	for (i = 0; i < n; ++i) D[i][i] = pow(D[i][i], -0.5);
}
double ** left_multip_of_diagonal_matrice(double ** D, double **M, int n) {
	int i, j;
	double **result;
	result = (double**)malloc(n * sizeof(double*));
	if (result == NULL) return NULL;
	for (i = 0; i < n; ++i) {
		result[i] = (double*)malloc (n * sizeof(double));
		if (result[i] == NULL) return NULL;
	}
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			result[i][j] = D[i][i] * M[i][j];
		}
	}
	return result;
}
double ** right_multip_of_diagonal_matrice(double ** M, double **D, int n) {
	int i, j;
	double **result;
	result = (double**)malloc(n * sizeof(double*));
	if (result == NULL) return NULL;
	for (i = 0; i < n; ++i) {
		result[i] = (double*)malloc(n * sizeof(double));
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
double ** calcNGL(double ** W, double ** D, int n) { 
	int i, j;
	double ** result, **tmp;
	raise_D_diagonal_in_minus_half(D, n);


	tmp = left_multip_of_diagonal_matrice(D, W, n);
	if (tmp == NULL) return NULL;
	result = right_multip_of_diagonal_matrice(tmp, D, n);
	if (result == NULL) return NULL;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			if (i == j) result[i][i] = 1 - result[i][i];
			else result[i][j] = -result[i][j];
		}
	}
	for (i = 0; i < n; ++i) free(tmp[i]);
	free(tmp);

	return result;
}