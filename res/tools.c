#ifndef TOOLS_C
#define TOOLS_C
#include "definitions.h"
#include "math.h"
double distance(double *v, double *u, int vectorDimension) {
	int i;
	double distance;
	distance = 0;
	for (i=0; i<vectorDimension; ++i) distance += pow((v[i] - u[i]), 2);
	return pow(distance, 0.5);
}
#endif