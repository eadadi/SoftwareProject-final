#ifndef TOOLS_C
#define TOOLS_C
#include <stdio.h>
static double distance(double *v, double *u, int vectorDimension) {
	int d = vectorDimension;
	double distance = 0;
	for (d; d > 0; --d) distance += pow((v[d - 1] - u[d - 1]), 2);
	return pow(distance, 0.5);
}

static void Transpoze(double ** V, int n) {
	int i, j;
	double tmpval;
	for (i = 0; i < n; ++i) {
		for (j = i + 1; j < n; ++j) {
			tmpval = V[i][j];
			V[i][j] = V[j][i];
			V[j][i] = tmpval;
		}
	}
}

static void printArr(double ** A, int n, int m) {
	int i, j;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < m; ++j) {
			printf("%lf", A[i][j]);
			j + 1 == m ? printf("\n") : printf(",");
		}
		if (m == 0) {
			printf("%lf", *A[i]);
			i + 1 == n ? printf("\n") : printf(",");
		}
	}
}
#endif