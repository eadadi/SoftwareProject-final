#include <math.h>
#include <stdlib.h>
#include<stdio.h>
#define EPSILON 0.001
#define JACOBI_MAX_ITERATIONS_NUMBER 100
double*** calcJacobi(double **, int);
double ** buildPij(double **, int, int, int, double, double);
int* getPivotIndexes(double **, int);
double* calcCandS(double **, int, int);
double ** MatriceMultiplication(double **, double **, int);
double ** calcAtag(double **, int, int, int, double, double);
int isDiagonalEnough(double **A, double **, int);

void printArr(double **A, int n) {
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j) {
			printf("%lf", A[i][j]); j + 1 == n ? printf("\n") : printf(",");
		}
}

double*** calcJacobi(double **A, int n) {
	double **Pij, **Atag, **tmp, **V, **Acopy;
	double ***result = (double***)calloc(2,sizeof(double**)); //Composed of eigenvectors and eigenvalues of A

	Acopy = (double**)calloc(n, sizeof(double*));
	for (int i = 0; i < n; ++i) {
		Acopy[i] = (double*)calloc(n, sizeof(double));
		for (int j = 0; j < n; ++j) Acopy[i][j] = A[i][j];
	}
	printf("A:\n");
	printArr(A,n);
	printf("Acopy:\n");
	printArr(Acopy,n);

	V = (double**)calloc(n, sizeof(double*));
	for (int i = 0; i < n; ++i) {
		V[i] = (double*)calloc(n, sizeof(double));
		V[i][i] = 1;
	}

	printf("V:\n");
	printArr(V,n);

	int iter = 0;
	while (iter < JACOBI_MAX_ITERATIONS_NUMBER) {
		int i, j;
		double c, s;
		int *pivotIndexes = getPivotIndexes(Acopy, n);
		i = pivotIndexes[0]; j = pivotIndexes[1];
		double *params = calcCandS(Acopy, i, j);
		c = params[0]; s = params[1];
		free(pivotIndexes); free(params);

		Pij = buildPij(Acopy, n, i, j, c, s);
		Atag = calcAtag(Acopy, n, i, j, c, s);

		
		

		if (isDiagonalEnough(Acopy, Atag, n)) {
			iter = JACOBI_MAX_ITERATIONS_NUMBER-1;
			printf("---iter:%d---\n", iter);
			printf("i=%d, j=%d, c=%lf, s=%lf\n", i, j, c, s);
			printf("Pij:\n");
			printArr(Pij, n);
			printf("Atag:\n");
			printArr(Atag, n);
		}
		else if (iter > JACOBI_MAX_ITERATIONS_NUMBER-3) {
			printf("---iter:%d---\n", iter);
			printf("i=%d, j=%d, c=%lf, s=%lf\n", i, j, c, s);
			printf("Pij:\n");
			printArr(Pij, n);
			printf("Atag:\n");
			printArr(Atag, n);
		}

		tmp = V;
		V = MatriceMultiplication(tmp, Pij, n);

		for (int k = 0; k < n; ++k) {
			free(tmp[k]); free(Pij[k]); free(Acopy[k]);
		}
		free(tmp); free(Pij); free(Acopy);

		Acopy = Atag;
		++iter;
	}

	double * eigenvalues = (double*)calloc(n, sizeof(double));
	for (int i = 0; i < n; ++i) eigenvalues[i] = Acopy[i][i];

	for (int i = 0; i < n; ++i) free(Acopy[i]);
	free(Acopy);

	result[0] = &eigenvalues; result[1] = V;
	return result;
}
double ** buildPij(double **M, int n, int i, int j, double c, double s) {
	double **Pij;;
	Pij = (double**)calloc(n, sizeof(double*));
	for (int k = 0; k < n; ++k) Pij[k] = (double*)calloc(n, sizeof(double));
	for (int k = 0; k < n; k++) Pij[k][k] = 1;
	Pij[i][i] = Pij[j][j] = c;
	Pij[i][j] = s;  Pij[j][i] = -s;

	return Pij;
}
int* getPivotIndexes(double ** M, int n) {
	int i, j;
	double max = 0, tmp;
	for (int k = 0; k < n; ++k)
		for (int l = 0; l < n; ++l) {
			if (l == k) continue;
			tmp = fabs(M[k][l]);
			if (tmp > max) {
				max = tmp;
				i = l; j = k;
			}
		}
	int * result = (double*)calloc(2, sizeof(double));

	if (i < j) {
		result[0] = i; result[1] = j;
	}
	else {
		result[0] = j; result[1] = i;
	} 
	return result;
}
double* calcCandS(double **M, int i, int j) {
	double theta = (M[j][j] - M[i][i]) / (2 * M[i][j]);
	int thetaSign;
	if (theta >= 0) thetaSign = 1;
	else thetaSign = 0;
	double t = thetaSign / (fabs(theta) + pow(pow(theta, 2) + 1, 0.5));
	double c = 1 / pow(pow(t, 2) + 1, 0.5);
	double s = t * c;
	double * result = (double*)calloc(2, sizeof(double));
	result[0] = c; result[1] = s;
	return result;
}
double ** MatriceMultiplication(double ** A, double **B, int n) {
	double ** C = (double **)calloc(n, sizeof(double*));
	for (int i = 0; i < n; ++i) {
		C[i] = (double*)calloc(n, sizeof(double));
		for (int j = 0; j < n; ++j) {
			for (int k = 0; k < n; ++k) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	};
	return C;
}

double ** calcAtag(double ** A, int n, int i, int j, double c, double s) {
	double **Atag;
	Atag = (double**)calloc(n, sizeof(double*));
	for (int i = 0; i < n; ++i) {
		Atag[i] = (double*)calloc(n, sizeof(double));
		for (int j = 0; j < n; ++j) Atag[i][j] = A[i][j];
	}
	for (int r = 0; r < n; ++r) {
		Atag[r][i] = Atag[i][r] =  c * A[r][i] - s * A[r][j];
		Atag[r][j] = Atag[j][r] = c * A[r][j] + s * A[r][i];
	}
	double s2 = pow(s, 2), c2 = pow(c, 2),
		aii = A[i][i], ajj = A[j][j], aij = A[i][j];
	Atag[i][i] = c2 * aii + s2 * ajj - 2 * s*c*aij;
	Atag[j][j] = s2 * aii + c2 * ajj + 2 * s*c*aij;
	Atag[i][j] = Atag[j][i] = 0;
	return Atag;
}

int isDiagonalEnough(double **A, double **B, int n) {
	double offA = 0, offB = 0;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (j == i) continue;
			offA += pow(A[i][j], 2);
			offB += pow(B[i][j], 2);
		}
	}
	if (fabs(offA - offB) <= EPSILON) return 1;
	else return 0;
}

