#include <stdio.h>
#include <stdlib.h>
#include<math.h>
/*
-----jacobi.c headers
*/
#ifndef JACOBI_CONSTANTS
#define JACOBI_CONSTANTS
#define EPSILON 0.001
#define JACOBI_MAX_ITERATIONS_NUMBER 100
#endif 


/*
-----matrice_max_heap.c headers
*/
#ifndef MATRICE_MAX_HEAP
#define MATRICE_MAX_HEAP
typedef struct matrice_max_heap {
	double * values;
	int **mat_to_values;
	int *values_to_mat;
}matrice_max_heap;
#define HEAP_MEM 2
#endif 

/*
-----value_vector_map_struct definition:
*/
#ifndef VALUE_VECTOR_STRUCT
#define VALUE_VECTOR_STRUCT
typedef struct value_vector_map {
	double **eigenvector;
	double *eigenvalue;
	int index;
}value_vector_map;
#endif

#ifndef C_ENUM_GOAL
#define C_ENUM_GOAL
enum goal { wam, ddg, lnorm, jacobi, spk };
#endif