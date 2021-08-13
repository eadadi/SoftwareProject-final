#ifndef MATRICE_MAX_HEAP
#define MATRICE_MAX_HEAP
typedef struct matrice_max_heap {
	double * values;
	int **mat_to_values;
	int *values_to_mat;
}matrice_max_heap;
#endif 