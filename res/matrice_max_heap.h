#ifndef JACOBI_CONSTANTS
#define JACOBI_CONSTANTS
#define EPSILON 0.001
#define JACOBI_MAX_ITERATIONS_NUMBER 100
#define HEAP_MEM 2
typedef struct matrice_max_heap {
	double * values;
	int **mat_to_values;
	int *values_to_mat;
}matrice_max_heap;
#endif 