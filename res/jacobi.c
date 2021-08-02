#include <stdlib.h>
#include "tools.c"

#ifndef JACOBI_CONSTANTSS
#define EPSILON 0.001
#define JACOBI_MAX_ITERATIONS_NUMBER 100
#endif 


static int* getPivotIndexes(double ** M, int n) {
	int i, j, k, l, *result;
	double max, tmp;
	result = (int*)calloc(2, sizeof(int));
	if (result == NULL) return NULL;
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
static double* calcCandS(double **M, int i, int j) {
	double theta, t, c, s, *result;
	int thetaSign;
	result = (double*)calloc(2, sizeof(double));
	if (result == NULL) return NULL;
	theta = (M[j][j] - M[i][i]) / (2 * M[i][j]);
	if (theta >= 0) thetaSign = 1;
	else thetaSign = -1;
	t = thetaSign / (fabs(theta) + pow(pow(theta, 2) + 1, 0.5));
	c = 1 / pow(pow(t, 2) + 1, 0.5);
	s = t * c;
	result[0] = c; result[1] = s;
	return result;
}
static double ** calcAtag(double ** A, int n, int i, int j, double c, double s, double *offAtag) {
	double c2, s2, aii, ajj, aij, ari, arj, **Atag;
	int r;
	Atag = A;
	c2 = pow(c, 2), s2 = pow(s, 2), aii = A[i][i], ajj = A[j][j], aij = A[i][j];
	for (r = 0; r < n; ++r) {
		ari = A[r][i];
		arj = A[r][j];
		*offAtag -= 2 * pow(ari, 2);
		*offAtag -= 2 * pow(arj, 2);
		Atag[r][i] = Atag[i][r] = c * ari - s * arj;
		Atag[r][j] = Atag[j][r] = c * arj + s * ari;
		*offAtag += 2 * pow(Atag[r][i], 2);
		*offAtag += 2 * pow(Atag[r][j], 2);
	}
	Atag[i][i] = c2 * aii + s2 * ajj - 2 * s*c*aij;
	Atag[j][j] = s2 * aii + c2 * ajj + 2 * s*c*aij;
	*offAtag -= 2 * pow(aij, 2);
	Atag[i][j] = Atag[j][i] = 0;
	return Atag;
}
static void update_V_by_Pij(double ** V, int n, int i, int j, double c, double s) {
	int r;
	double Vri, Vrj;
	for (r = 0; r < n; ++r) {
		Vri = V[r][i];
		Vrj = V[r][j];
		V[r][i] = c * Vri - s * Vrj;
		V[r][j] = s * Vri + c * Vrj;
	}
}
static double calcOffA(double ** A, int n) {
	int i, j;
	double offA = 0;
	for (i = 0; i < n; ++i) {
		for (j = i+1; j < n; ++j) {
			offA += 2 * pow(A[i][j], 2);
		}
	}
	return offA;
}
static double*** calcJacobi(double **A, int n) {
	double **Atag, **V, **Acopy, ***result;
	int iter, i, j, *pivotIndexes, flag;
	double c, s, *params, *eigenvalues, **eigenvaluesHolder;
	double offAcopy, offAtag;
	eigenvalues = (double*)calloc(n, sizeof(double));
	if (eigenvalues == NULL) return NULL;
	eigenvaluesHolder = (double**)calloc(1, sizeof(double*));
	if (eigenvaluesHolder == NULL) return NULL;
	result = (double***)calloc(2, sizeof(double**));
	if (result == NULL) return NULL;

	iter = 0;
	Acopy = (double**)calloc(n, sizeof(double*));
	V = (double**)calloc(n, sizeof(double*));
	if (Acopy == NULL || V == NULL) return NULL;
	for (i = 0; i < n; ++i) {
		Acopy[i] = (double*)calloc(n, sizeof(double));
		if (Acopy[i] == NULL) return NULL;
		V[i] = (double*)calloc(n, sizeof(double));
		if (V[i] == NULL) return NULL;
	}
	for (i = 0; i < n; ++i) {
		V[i][i] = 1;
		for (j = 0; j < n; ++j) {
			Acopy[i][j] = A[i][j];
		}
	}

	offAcopy = calcOffA(Acopy, n);
	flag = 1;
	while (flag && iter++ < JACOBI_MAX_ITERATIONS_NUMBER) {
		pivotIndexes = getPivotIndexes(Acopy, n);
		if (pivotIndexes == NULL) return NULL;
		i = pivotIndexes[0]; j = pivotIndexes[1];
		params = calcCandS(Acopy, i, j);
		if (params == NULL) return NULL;
		c = params[0]; s = params[1];
		free(pivotIndexes); free(params);

		offAtag = offAcopy;

		/* calcAtag returns the Acopy ptr (it updates Acopy) */
		Atag = calcAtag(Acopy, n, i, j, c, s, &offAtag); 

		if (Atag == NULL) return NULL;
		if (offAcopy - offAtag <= EPSILON) flag = 0;
		update_V_by_Pij(V, n, i, j, c, s);
		Acopy = Atag;
		offAcopy = offAtag;
	}

	for (i = 0; i < n; ++i) {
		eigenvalues[i] = Acopy[i][i];
		free(Acopy[i]);
	}free(Acopy);
	
	Transpoze(V, n);

	eigenvaluesHolder[0] = eigenvalues;
	result[0] = eigenvaluesHolder; result[1] = V;
	return result;
}