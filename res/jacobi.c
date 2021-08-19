#include <stdlib.h>
#include "tools.c"
#include <float.h>
#include "../spkmeans.h"
#include "matrice_max_heap.h"

static int parent_loc(int i) {
	int l;
	l = i / 2;
	return l;
}
static double parent(matrice_max_heap *h, int i) {
	int l;
	l = parent_loc(i);
	return h->values[l];
}
static int left_loc(int i) {
	int l;
	l = 2 * i;
	return l;
}
static double left(matrice_max_heap *h, int i) {
	int l;
	l = left_loc(i);
	return h->values[l];
}
static int right_loc(int i) {
	int l;
	l = 2 * i + 1;
	return l;
}
static double right(matrice_max_heap *h, int i) {
	int l;
	l = right_loc(i);
	return h->values[l];
}
static double key(matrice_max_heap *h, int i) {
	return h->values[i];
}
static void set_key(matrice_max_heap *h, int i, double k) {
	h->values[i] = k;
}
static void set_indexes(matrice_max_heap *h, int index, int i, int j) {
	h->values_to_mat[index * HEAP_MEM] = i;
	h->values_to_mat[index * HEAP_MEM + 1] = j;
}
static int heap_len(matrice_max_heap *h) {
	return h->values_to_mat[0];
}
static void interswitch(matrice_max_heap *h, int i, int j) {
	int j_ind0, j_ind1, i_ind0, i_ind1, tmp;
	double j_k, i_k;
	j_k = key(h, j);
	j_ind0 = h->values_to_mat[HEAP_MEM * j];
	j_ind1 = h->values_to_mat[HEAP_MEM * j + 1];

	i_k = key(h, i);
	i_ind0 = h->values_to_mat[HEAP_MEM * i];
	i_ind1 = h->values_to_mat[HEAP_MEM * i + 1];

	tmp = h->mat_to_values[j_ind0][j_ind1];
	h->mat_to_values[j_ind0][j_ind1] = h->mat_to_values[i_ind0][i_ind1];
	h->mat_to_values[i_ind0][i_ind1] = tmp;

	set_key(h, j, i_k);
	set_indexes(h, j, i_ind0, i_ind1);
	set_key(h, i, j_k);
	set_indexes(h, i, j_ind0, j_ind1);
}
static void heapify_down(matrice_max_heap *h, int i) {
	int l, r, max;
	l = left_loc(i);
	r = right_loc(i);
	max = i;
	if (l<heap_len(h) && left(h,i)>key(h, max))
		max = l;
	if (r<heap_len(h) && right(h,i)>key(h, max))
		max = r;
	if (max > i) {
		interswitch(h, i, max);
		heapify_down(h, max);
	}
}
static void heapify_up(matrice_max_heap *h, int i) {
	while (i > 1 && key(h, i) > parent(h, i)) {
		interswitch(h, i, parent_loc(i));
		i = parent_loc(i);
	}
}
static void heapify(matrice_max_heap *h) {
	int i, l;
	l = heap_len(h);
	for (i = l-1; i >= 1; --i) {
		heapify_down(h, i);
	}
}

/*[i, j] where max value obtained in a_ij*/
static int* heap_max(matrice_max_heap *h) { 
	int * result, i, j;
	result = (int*)malloc(2 * sizeof(int));
	
	i = h->values_to_mat[1 * HEAP_MEM];
	j = h->values_to_mat[1 * HEAP_MEM + 1];
	result[0] = i;
	result[1] = j;

	return result;
}

static void update_key(matrice_max_heap *h, int _ind, double v) {
	if (fabs(v) > key(h, _ind)) {
		set_key(h, _ind, fabs(v));
		heapify_up(h, _ind);
	}
	else if (fabs(v) < key(h, _ind)) {
		set_key(h, _ind, fabs(v));
		heapify_down(h, _ind);
	}
}
/*return -1 on failure*/
static int init_max_heap_from_matrice(matrice_max_heap *h, double **mat,
	int **heap_locations_map, int n) {
	int i, j, k, m;
	m = pow(n, 2) - n + 1;
	h->values = (double*)malloc(m * sizeof(double));
	if (h->values == NULL) return -1;
	h->values_to_mat = (int*)malloc(m * HEAP_MEM * sizeof(int));
	if (h->values_to_mat == NULL) return -1;

	h->mat_to_values = heap_locations_map;
	h->values[0] = -1; 
	h->values_to_mat[0] = m;
	h->values_to_mat[1] = -1;
	k = 1;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			if (i == j)continue;
			h->values[k] = fabs(mat[i][j]);
			h->values_to_mat[k * HEAP_MEM] = i;
			h->values_to_mat[k * HEAP_MEM + 1] = j;
			heap_locations_map[i][j] = k;
			++k;
		}
	}
	heapify(h);

	
	return 0;
}

static void free_heap(matrice_max_heap *h, int mat_dim) {
	int i;
	free(h->values);
	free(h->values_to_mat);
	for (i = 0; i < mat_dim; ++i) free(h->mat_to_values[i]);
	free(h->mat_to_values);
}

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
static void calcCandS(double **M, int i, int j, double result[2]) {
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
static double ** calcAtag(double ** A, int n, int i, int j, 
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
		for (j = 0; j < n; ++j) {
			if (i == j) continue;
			offA += pow(A[i][j], 2);
		}
	}
	return offA;
}
static double*** calcJacobi(double **A, int n) {
	double **Atag, **V, ***result;
	int iter, i, j, pivotIndexes[2], flag;
	double c, s, params[2], *eigenvalues, **eigenvaluesHolder;
	double offA, offAtag;
	matrice_max_heap h;/**/
	int *max;
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
		 getPivotIndexes(A, n, pivotIndexes);
		i = pivotIndexes[0]; j = pivotIndexes[1];
		calcCandS(A, i, j, params);
		c = params[0]; s = params[1];
		max = heap_max(&h);

		free(max);

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
