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
static double ** raise_D_diagonal_in_minus_half(double **D, int n);
static double ** left_multip_of_diagonal_matrice(double ** D, double **M, int n);
static double ** right_multip_of_diagonal_matrice(double ** M, double **D, int n);
static double ** calcNGL(double ** W, double ** D, int n);

/*
-----jacobi.c headers
*/
#ifndef JACOBI_CONSTANTSS
#define EPSILON 0.001
#define JACOBI_MAX_ITERATIONS_NUMBER 100
#endif // !JACOBI_CONSTANTSS
static int* getPivotIndexes(double ** M, int n);
static double* calcCandS(double **M, int i, int j);
static double ** calcAtag(double ** A, int n, int i, int j, double c, double s, double *offAtag);
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



#endif