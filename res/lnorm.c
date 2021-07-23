#include <math.h>
#include <stdlib.h>
double ** raise_D_diagonal_in_minus_half(double **D, int n) {
	double ** raised = (double**)calloc(n, sizeof(double*));
	for (int i = 0; i < n; ++i) raised[i] = (double*)calloc(n, sizeof(double));
	for (int i = 0; i < n; ++i) raised[i][i] = pow(D[i][i], -0.5);
	return raised;
}
double ** left_multip_of_diagonal_matrice(double ** D, double **M, int n) {
	double ** result = (double**)calloc(n, sizeof(double*));
	for (int i = 0; i < n; ++i) result[i] = (double*)calloc(n, sizeof(double));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			result[i][j] = D[i][i] * M[i][j];
		}
	}
	return result;
}
double ** right_multip_of_diagonal_matrice(double ** M, double **D, int n) {
	double ** result = (double**)calloc(n, sizeof(double*));
	for (int i = 0; i < n; ++i) result[i] = (double*)calloc(n, sizeof(double));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			result[i][j] = M[i][j] * D[j][j];
		}
	}
	return result;
}
double ** calcNGL(double ** W, double ** D, int n) { // Normalized Graph Laplacian
	double ** result;
	double ** Dtag = raise_D_diagonal_in_minus_half(D, n);
	double ** tmp = left_multip_of_diagonal_matrice(Dtag, W, n);
	result = right_multip_of_diagonal_matrice(tmp, Dtag, n);
	for (int i = 0; i < n; ++i) 
		for (int j = 0; j < n; ++j)
			if (i == j) result[i][i] = 1 - result[i][i];
			else result[i][j] = -result[i][j];
	
	for (int i = 0; i < n; ++i) {
		free(Dtag[i]); free(tmp[i]); \
	}
	free(Dtag); free(tmp);
	
	return result;
}