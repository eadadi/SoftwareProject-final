#include <stdlib.h>
#include "tools.c"
static void add_u_to_v(double *u, double *v, int dim) {
	int i;
	for (i = 0; i < dim; ++i) {
		v[i] += u[i];
	}
}

static int updateCentroids(double **datapoints, double **centroids,
	int datapoints_amount, int clusters_amount, int *data_to_centroids_map, int datapoint_length) {
	double **new_centroids;
	int *amounts, i, j, flag;

	amounts = (int*)calloc(clusters_amount, sizeof(int));
	new_centroids = (double**)calloc(clusters_amount, sizeof(double*));
	if (amounts == NULL || new_centroids == NULL) return -1;
	for (i = 0; i < clusters_amount; ++i) {
		new_centroids[i] = (double*)calloc(datapoint_length, sizeof(double));
		if (new_centroids[i] == NULL) return -1;
	}
	for (i = 0; i < datapoints_amount; ++i) {
		j = data_to_centroids_map[i];
		add_u_to_v(datapoints[i], new_centroids[j], datapoint_length);
		amounts[j]++;
	}
	flag = 0;
	for (i = 0; i < clusters_amount; ++i) {
		for (j = 0; j < datapoint_length; ++j) {
			new_centroids[i][j] /= amounts[i];
			if (new_centroids[i][j] != centroids[i][j]) {
				centroids[i][j] = new_centroids[i][j];
				flag = 1;
			}
		}
	}
	//free new_centroids
	for (i = 0; i < clusters_amount; ++i) free(new_centroids[i]);
	free(new_centroids);
	//free amounts
	free(amounts);
	return flag;
}

static int getClosestCluster(double *candidate, double**centroids, int clusters_amount,
	int datapoint_length) {
	int i, index;
	double min, tmp;
	min = distance(candidate, centroids[0], datapoint_length);
	index = 0;
	for (i = 1; i < clusters_amount; ++i) {
		tmp = distance(candidate, centroids[i], datapoint_length);
		if (tmp < min) {
			min = tmp;
			index = i;
		}
	}
	return index;
}

static int* k_mean(double **datapoints, double **centroids, int datapoints_amount,
	int clusters_amount, int datapoint_length, int max_iter) {
	int *data_to_centroids_map, iter, i, closest_cluster, flag;
	flag = 1; iter = 0;
	data_to_centroids_map = (int*)malloc(datapoints_amount * sizeof(int));
	if (data_to_centroids_map == NULL) return NULL;
	while (flag == 1 && iter++ < max_iter) {
		for (i = 0; i < datapoints_amount; ++i) {
			closest_cluster = getClosestCluster(datapoints[i], centroids, clusters_amount, datapoint_length);
			data_to_centroids_map[i] = closest_cluster;
		}
		flag = updateCentroids(datapoints, centroids,
			datapoints_amount, clusters_amount, data_to_centroids_map, datapoint_length);
		if (flag == -1) return NULL;
	}
	return data_to_centroids_map;
}