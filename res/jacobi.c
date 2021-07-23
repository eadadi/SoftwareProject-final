#include <math.h>
#include <stdlib.h>
//#include<stdio.h>
#define EPSILON 0.001
#define JACOBI_MAX_ITERATIONS_NUMBER 100

double*** calcJacobi(double **, int);
double ** buildPij(double **, int, int, int, double, double);
int* getPivotIndexes(double **, int);
double* calcCandS(double **, int, int);
double ** MatriceMultiplication(double **, double **, int);
double ** calcAtag(double **, int, int, int, double, double);
int isDiagonalEnough(double **, double **, int);

double*** calcJacobi(double **A, int n) {
	double **Pij, **Atag, **tmp, **V, **Acopy, ***result;
	int iter = 0, i, j, *pivotIndexes;
	double c, s, *params, *eigenvalues, **eigenvaluesHolder;

	Acopy = (double**)calloc(n, sizeof(double*));
	V = (double**)calloc(n, sizeof(double*));
	for (int i = 0; i < n; ++i) {
		Acopy[i] = (double*)calloc(n, sizeof(double));
		for (int j = 0; j < n; ++j) Acopy[i][j] = A[i][j];

		V[i] = (double*)calloc(n, sizeof(double));
		V[i][i] = 1;
	}
	while (iter++ < JACOBI_MAX_ITERATIONS_NUMBER) {
		pivotIndexes = getPivotIndexes(Acopy, n);
		i = pivotIndexes[0]; j = pivotIndexes[1];
		params = calcCandS(Acopy, i, j);
		c = params[0]; s = params[1];
		free(pivotIndexes); free(params);

		Pij = buildPij(Acopy, n, i, j, c, s);
		Atag = calcAtag(Acopy, n, i, j, c, s);

		tmp = V; V = MatriceMultiplication(tmp, Pij, n);

		for (int k = 0; k < n; ++k) {
			free(tmp[k]); free(Pij[k]); free(Acopy[k]);
		} free(tmp); free(Pij); free(Acopy);
		Acopy = Atag;
	}

	eigenvalues = (double*)calloc(n, sizeof(double));
	for (int i = 0; i < n; ++i) { 
		eigenvalues[i] = Acopy[i][i]; 
		free(Acopy[i]);
	}free(Acopy);
	eigenvaluesHolder = (double**)calloc(1, sizeof(double*));
	eigenvaluesHolder[0] = eigenvalues;

	result = (double***)calloc(2, sizeof(double**));
	result[0] = eigenvaluesHolder; result[1] = V;
	return result;
}
double ** buildPij(double **M, int n, int i, int j, double c, double s) {
	double **Pij;;
	Pij = (double**)calloc(n, sizeof(double*));
	for (int k = 0; k < n; ++k) {
		Pij[k] = (double*)calloc(n, sizeof(double));
		Pij[k][k] = 1;
	}
	Pij[i][i] = Pij[j][j] = c;
	Pij[i][j] = s;  Pij[j][i] = -s;
	return Pij;
}
int* getPivotIndexes(double ** M, int n) {
	int i, j, k, l, *result;
	double max, tmp;
	result = (int*)calloc(2, sizeof(int));
	max = 0;
	for (k = 0; k < n; ++k)
		for (l = 0; l < n; ++l) {
			if (l == k) continue;
			tmp = fabs(M[k][l]);
			if (tmp > max) {
				max = tmp;
				i = l; j = k;
			}
		}
	result[0] = i; result[1] = j;
	return result;
}
double* calcCandS(double **M, int i, int j) {
	double theta, t, c, s, *result;
	int thetaSign;
	result = (double*)calloc(2, sizeof(double));
	theta = (M[j][j] - M[i][i]) / (2 * M[i][j]);
	if (theta >= 0) thetaSign = 1;
	else thetaSign = -1;
	t = thetaSign / (fabs(theta) + pow(pow(theta, 2) + 1, 0.5));
	c = 1 / pow(pow(t, 2) + 1, 0.5);
	s = t * c;
	result[0] = c; result[1] = s;
	return result;
}
double ** MatriceMultiplication(double ** A, double **B, int n) {
	double ** C;
	int i, j, k;
	C = (double **)calloc(n, sizeof(double*));
	for (i = 0; i < n; ++i) {
		C[i] = (double*)calloc(n, sizeof(double));
		for (j = 0; j < n; ++j) {
			for (k = 0; k < n; ++k) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	};
	return C;
}
double ** calcAtag(double ** A, int n, int i, int j, double c, double s) {
	double c2, s2, aii, ajj, aij, **Atag;
	int k, l, r;
	Atag = (double**)calloc(n, sizeof(double*));
	for (k = 0; k < n; ++k) {
		Atag[k] = (double*)calloc(n, sizeof(double));
		for (l = 0; l < n; ++l) Atag[k][l] = A[k][l];
	}
	for (r = 0; r < n; ++r) {
		Atag[r][i] = Atag[i][r] =  c * A[r][i] - s * A[r][j];
		Atag[r][j] = Atag[j][r] = c * A[r][j] + s * A[r][i];
	}
	c2 = pow(c, 2), s2 = pow(s, 2), aii = A[i][i], ajj = A[j][j], aij = A[i][j];
	Atag[i][i] = c2 * aii + s2 * ajj - 2 * s*c*aij;
	Atag[j][j] = s2 * aii + c2 * ajj + 2 * s*c*aij;
	Atag[i][j] = Atag[j][i] = 0;
	return Atag;
}
int isDiagonalEnough(double **A, double **B, int n) {
	double offA, offB;
	int i, j;
	offA = offB = 0;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			if (j == i) continue;
			offA += pow(A[i][j], 2);
			offB += pow(B[i][j], 2);
		}
	}
	if (fabs(offA - offB) <= EPSILON) return 1;
	else return 0;
}

