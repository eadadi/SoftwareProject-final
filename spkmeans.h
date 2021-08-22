#ifndef SPKMEANS_H
#define SPKMEANS_H

#include "res/definitions.h"
/*
-------------spkmeansmodule.c sources:
*/

/*
-----tools.c headers
*/
double distance(double *v, double *u, int vectorDimension);



/*
-----wam.c, ddg.c headers
*/
double ** calcWAM(double ** vectors, int n, int vectorDimension);
double **calcDDG(double ** M, int n);

/*
-----lnorm.c headers
*/
void raise_D_diagonal_in_minus_half(double **D, int n);
double ** left_multip_of_diagonal_matrice(double ** D, double **M, int n);
double ** right_multip_of_diagonal_matrice(double ** M, double **D, int n);
double ** calcNGL(double ** W, double ** D, int n);

/*
-----matrice_max_heap.c headers
*/

int parent_loc(int);
double parent(matrice_max_heap*, int);
int left_loc(int i);
double left(matrice_max_heap*, int);
int right_loc(int i);
double right(matrice_max_heap*, int);
double key(matrice_max_heap*, int);
void set_key(matrice_max_heap*, int, double);
void set_indexes(matrice_max_heap*, int, int, int);
int heap_len(matrice_max_heap*);
void interswitch(matrice_max_heap*, int, int);
void heapify_down(matrice_max_heap*, int);
void heapify_up(matrice_max_heap*, int);
void heapify(matrice_max_heap*);
void heap_max(matrice_max_heap*, int[2]);
void update_key(matrice_max_heap*, int, double);
int init_max_heap_from_matrice(matrice_max_heap*, double **, int **, int);
void free_heap(matrice_max_heap*, int);

/*
-----jacobi.c headers
*/

/*static void getPivotIndexes(double ** M, int n, int[2]);*/
void calcCandS(double **M, int i, int j, double[2]);
double ** calcAtag(double ** A, int n, int i, int j, double c, 
	double s, double *offAtag, matrice_max_heap *h);
void update_V_by_Pij(double ** V, int n, int i, int j, double c, double s);
double calcOffA(double ** A, int n);
void Transpoze(double ** V, int n);
double*** calcJacobi(double **A, int n);

/*
-----eigenpap.c headers
*/
int cmpfunc(const void *a, const void *b);
void sortMap(value_vector_map *map, int n);
int determineK(value_vector_map *map, int n);
value_vector_map* setMap(double ** eigenvectors, double * eigenvalues, int n);

/*
-----fit.c headers
*/
void add_u_to_v(double *u, double *v, int dim);
int updateCentroids(double **datapoints, double **centroids,
	int datapoints_amount, int clusters_amount, int *data_to_centroids_map, int datapoint_length);
int getClosestCluster(double *candidate, double**centroids, int clusters_amount,
	int datapoint_length);
int* k_mean(double **datapoints, double **centroids, int datapoints_amount,
	int clusters_amount, int datapoint_length, int max_iter);
/*
-----spk.c headers
*/
double ** calcNormalaizedRnk(value_vector_map *map, int n, int clustersNumber);
void printArr(double ** A, int n, int m);


/*
-------------spkmeans.c  headers:
*/
double **extractCentroids(double **data, int k);
double ** calcFinalCentroids(double **data, int *map, int k, int data_length, int data_dim);

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