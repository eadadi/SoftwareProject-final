#ifndef SPKMEANS_H
#define SPKMEANS_H

/*
-------------spkmeansmodule.c sources:
*/

/*
-----tools.c headers
*/
static double distance(double *v, double *u, int vectorDimension);
static void Transpoze(double ** V, int n);
static void printArr(double ** A, int n, int m);

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

/*
-----wam.c, ddg.c headers
*/
static double ** calcWAM(double ** vectors, int n, int vectorDimension);
static double **calcDDG(double ** M, int n);

/*
-----lnorm.c headers
*/
static void raise_D_diagonal_in_minus_half(double **D, int n);
static double ** left_multip_of_diagonal_matrice(double ** D, double **M, int n);
static double ** right_multip_of_diagonal_matrice(double ** M, double **D, int n);
static double ** calcNGL(double ** W, double ** D, int n);

/*
-----jacobi.c headers
*/
#ifndef JACOBI_CONSTANTS
#define JACOBI_CONSTANTS
#define EPSILON 0.0001
#define JACOBI_MAX_ITERATIONS_NUMBER 10000
#define HEAP_MEM 2
#endif 
#ifndef MATRICE_MAX_HEAP
#define MATRICE_MAX_HEAP
typedef struct matrice_max_heap {
	double * values;
	int **mat_to_values;
	int *values_to_mat;
}matrice_max_heap;
#endif 

static void getPivotIndexes(double ** M, int n, int[2]);
static void calcCandS(double **M, int i, int j, double[2]);
static double ** calcAtag(double ** A, int n, int i, int j, double c, 
	double s, double *offAtag, matrice_max_heap *h);
static void update_V_by_Pij(double ** V, int n, int i, int j, double c, double s);
static double calcOffA(double ** A, int n);
static double*** calcJacobi(double **A, int n);

/*
-----eigenpap.c headers
*/
static int cmpfunc(const void *a, const void *b);
static void sortMap(value_vector_map *map, int n);
static int determineK(value_vector_map *map, int n);
static value_vector_map* setMap(double ** eigenvectors, double * eigenvalues, int n);

/*
-----fit.c headers
*/
static void add_u_to_v(double *u, double *v, int dim);
static int updateCentroids(double **datapoints, double **centroids,
	int datapoints_amount, int clusters_amount, int *data_to_centroids_map, int datapoint_length);
static int getClosestCluster(double *candidate, double**centroids, int clusters_amount,
	int datapoint_length);
static int* k_mean(double **datapoints, double **centroids, int datapoints_amount,
	int clusters_amount, int datapoint_length, int max_iter);
/*
-----spk.c headers
*/
static double ** calcNormalaizedRnk(value_vector_map *map, int n, int clustersNumber);


/*
-------------spkmeans.c  headers:
*/
double **extractCentroids(double **data, int k);
double ** calcFinalCentroids(double **data, int *map, int k, int data_length, int data_dim);

#ifndef C_ENUM_GOAL
#define C_ENUM_GOAL
enum goal { wam, ddg, lnorm, jacobi, spk };
#endif

char* __goal(enum goal e);
int cmpstr(char * a, char *b);
int lenstr(char * a);
void copyAtoB(char *A, char *B);
int freadline(char **a, FILE *f);
int countCommas(char *a);
double *commaSplit(char *line, int features);
int getdata(FILE *f, double ***data, int *features_num);
int getgoal(enum goal* e, char *candidate);
int input(int argc, char*argv[], int *k, enum goal* e, double *** data, int *data_len, int *features_number);

#endif