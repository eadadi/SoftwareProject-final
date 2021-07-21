#include <math.h>
#include <stdlib.h>
double distance(double *v, double *u, int vectorDimension) {
	int d = vectorDimension;
	double distance = 0;
	for (d; d > 0; --d) distance += pow((v[d - 1] - u[d - 1]), 2);
	return pow(distance, 0.5);
}
double ** calcWAM (double ** vectors, int n, int vectorDimension){
    int d = vectorDimension;
    double ** Matrice = (double**)calloc(n,sizeof(double*));
	for (int i = 0; i < n; ++i){
		Matrice[i] = (double*)calloc(n,sizeof(double));
	}
	for (int i = 0; i < n; ++i) {
		for (int j = i; j < n; ++j) {
			if (j == i) Matrice[i][j] = 0;
			else Matrice[i][j] = Matrice[j][i] = distance(vectors[i], vectors[j], d);
		}
	}
    return Matrice;
}

