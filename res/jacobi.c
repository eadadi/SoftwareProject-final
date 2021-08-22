#include "../spkmeans.h"
#include "matrice_max_heap.c"
/*
static void getPivotIndexes(double ** M, int n, int result[2]) {
	int i, j, k, l;
	double max, tmp;
	max = fabs(M[0][1]);
	i = 0; j = 1;
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
}
*/
void calcCandS(double **M, int i, int j, double result[2]) {
	double theta, t, c, s;
	int thetaSign;
	theta = (M[j][j] - M[i][i]) / (2 * M[i][j]);
	if (theta >= 0) thetaSign = 1;
	else thetaSign = -1;
	t = thetaSign / (fabs(theta) + pow(pow(theta, 2) + 1, 0.5));
	c = 1 / pow(pow(t, 2) + 1, 0.5);
	s = t * c;
	result[0] = c; result[1] = s;
}


double ** calcAtag(double ** A, int n, int i, int j, 
	double c, double s, double *offAtag, matrice_max_heap *h) {
	double c2, s2, aii, ajj, aij, ari, arj, **Atag, v1, v2;
	int r;
	Atag = A;
	c2 = pow(c, 2), s2 = pow(s, 2), aii = A[i][i], ajj = A[j][j], aij = A[i][j];
	for (r = 0; r < n; ++r) {
		if (r == i || r == j) continue;
		ari = A[r][i];
		arj = A[r][j];
		*offAtag -= 2 * pow(ari, 2);
		*offAtag -= 2 * pow(arj, 2);
		v1 = c * ari - s * arj;
		v2 = c * arj + s * ari;
		Atag[r][i] = Atag[i][r] = v1;
		Atag[r][j] = Atag[j][r] = v2;
		/*heapifing*/
		update_key(h, h->mat_to_values[i][r], v1);
		update_key(h, h->mat_to_values[r][i], v1);
		update_key(h, h->mat_to_values[j][r], v2);
		update_key(h, h->mat_to_values[r][j], v2);
		
		*offAtag += 2 * pow(v1, 2);
		*offAtag += 2 * pow(v2, 2);
	}
	Atag[i][i] = c2 * aii + s2 * ajj - 2 * s*c*aij;
	Atag[j][j] = s2 * aii + c2 * ajj + 2 * s*c*aij;
	*offAtag -= 2 * pow(aij, 2);
	Atag[i][j] = Atag[j][i] = 0;
	
	/*heapifing*/
	update_key(h, h->mat_to_values[i][j], 0);
	update_key(h, h->mat_to_values[j][i], 0);

	return Atag;
}


void update_V_by_Pij(double ** V, int n, int i, int j, double c, double s) {
	int r;
	double Vri, Vrj;
	for (r = 0; r < n; ++r) {
		Vri = V[r][i];
		Vrj = V[r][j];
		V[r][i] = c * Vri - s * Vrj;
		V[r][j] = s * Vri + c * Vrj;
	}
}


double calcOffA(double ** A, int n) {
	int i, j;
	double offA = 0;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			if (i == j) continue;
			offA += pow(A[i][j], 2);
		}
	}
	return offA;
}

void Transpoze(double ** V, int n) {
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

double*** calcJacobi(double **A, int n) {
	double **Atag, **V, ***result;
	int iter, i, j, flag;
	/*int pivotIndexes[2]; was for naive method of finding pivot */
	double c, s, params[2], *eigenvalues, **eigenvaluesHolder;
	double offA, offAtag;
	matrice_max_heap h;/* heap is now used to find pivot indexes*/
	int pivotIndexes[2];
	int **heap_locations;

	heap_locations = (int**)malloc(n * sizeof(int*));
	for (i = 0; i < n; ++i) heap_locations[i] = (int*)malloc(n * sizeof(int));


	eigenvalues = (double*)calloc(n, sizeof(double));
	if (eigenvalues == NULL) return NULL;
	eigenvaluesHolder = (double**)calloc(1, sizeof(double*));
	if (eigenvaluesHolder == NULL) return NULL;
	result = (double***)calloc(2, sizeof(double**));
	if (result == NULL) return NULL;

	iter = 0;
	V = (double**)calloc(n, sizeof(double*));
	if (V == NULL) return NULL;
	for (i = 0; i < n; ++i) {
		V[i] = (double*)calloc(n, sizeof(double));
		if (V[i] == NULL) return NULL;
	}
	for (i = 0; i < n; ++i) {
		V[i][i] = 1;
	}
	init_max_heap_from_matrice(&h, A,heap_locations, n);
	offA = calcOffA(A, n);
	flag = 1;
	
	while (flag==1 && iter++ < JACOBI_MAX_ITERATIONS_NUMBER) {
		/*getPivotIndexes(A, n, pivotIndexes);
		i = pivotIndexes[0]; j = pivotIndexes[1];*/
		heap_max(&h,pivotIndexes);
		i = pivotIndexes[0]; j = pivotIndexes[1];
		calcCandS(A, i, j, params);
		c = params[0]; s = params[1];
		

		offAtag = offA;

		/* calcAtag returns the A ptr (it updates A) */
		Atag = calcAtag(A, n, i, j, c, s, &offAtag, &h); 

		if (offA - offAtag <= EPSILON) flag = 0;

		update_V_by_Pij(V, n, i, j, c, s);

		A = Atag;


		offA = offAtag;
	}
	free_heap(&h, n);

	for (i = 0; i < n; ++i) {
		eigenvalues[i] = A[i][i];
	}
	Transpoze(V, n);
	eigenvaluesHolder[0] = eigenvalues;
	result[0] = eigenvaluesHolder; result[1] = V;

	return result;
}