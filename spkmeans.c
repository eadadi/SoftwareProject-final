#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "spkmeans.h"
#include "res/value_vector_map_struct.h"
#include "res/tools.c"
#include "res/wam.c"
#include "res/ddg.c"
#include "res/lnorm.c"
#include "res/jacobi.c"
#include "res/spk.c"
#include "res/eigenpap.c"
#include "res/fit.c"

#define NORMAL_LENGTH 10
#define MAX_NUM_OF_POINTS 1000

enum goal{wam, ddg, lnorm, jacobi, spk };
char* __goal(enum goal e) {
	switch (e)
	{
	case wam:
		return "wam";
	case ddg:
		return "ddg";
	case lnorm:
		return "lnorm";
	case jacobi:
		return "jacobi";
	case spk:
		return "spk";
	default:
		return "\0";
	}
}
int cmpstr(char * a, char *b) {
	int i;
	i = 0;
	if (a == NULL && b == NULL) return 1;
	if (a == NULL || b == NULL) return 0;
	while (a[i] != '\0' && b[i] != '\0') {
		if (a[i] != b[i]) return 0;
		++i;
	}
	if (a[i] == b[i]) return 1;
	return 0;
}
int lenstr(char * a) {
	int i;
	while (a[i++] != '\0');
	return i - 1;
}
void copyAtoB(char *A, char *B) {
	int k = 0;
	while ((B[k++] = A[k]) != '\0');
}
/*freadline returns -1 on failure*/
int freadline(char **a, FILE *f) { 
	char * tmp, ch, *backup;
	int i, curr_len = NORMAL_LENGTH;
	tmp = (char*)malloc(curr_len * sizeof(char));
	if (tmp == NULL) return -1;
	i = 0;
	while ((ch = fgetc(f)) != '\n' && ch != EOF) {
		tmp[i] = ch;
		if (i + 2 == curr_len) {
			tmp[i + 1] = '\0';
			curr_len *= 2;
			backup = (char*)malloc(curr_len * sizeof(char));
			if (backup == NULL) return -1;
			copyAtoB(tmp, backup);
			free(tmp);
			tmp = backup;
		}
		++i;
	}
	tmp[i] = '\0';
	*a = tmp;
	return i;
}
int countCommas(char *a) {
	int i,cnt;
	i = cnt = 0;
	while (a[i] != '\0') if (a[i++] == ',') ++cnt;
	return cnt;
}
/*commaSplit returns NULL on failure*/
double *commaSplit(char *line, int features) { 
	double * vector;
	char *ptr;
	int i;
	i = 0;
	vector = (double*)malloc(features * sizeof(double));
	if (vector == NULL) return NULL;
	while (i < features) {
		ptr = line + 1;
		while (*ptr != ',' && *ptr != '\0') ++ptr;
		*ptr = '\0';
		vector[i++] = atof(line);
		line = ptr + 1;
	}
	return vector;
}
/*getdata returns -1 on failure*/
int getdata(FILE *f,double ***data,int *features_num) { 
	double *vector;
	char *line;
	int i, len;
	*data = (double**)malloc(MAX_NUM_OF_POINTS * sizeof(double*));
	len = freadline(&line, f);
	if (len == -1) return -1;
	*features_num = countCommas(line) + 1;
	i = 0;
	while (line[0]!='\0') {
		vector = commaSplit(line, *features_num);
		if (vector == NULL) return -1;
		(*data)[i] = vector;
		free(line); line = NULL;
		len = freadline(&line, f);
		++i;
	}
	free(line);
	return i;
}
/* returns 0 on failure */
int getgoal(enum goal* e, char *candidate) {
	if (cmpstr(candidate, "wam")) *e = wam;
	else if (cmpstr(candidate, "ddg")) *e = ddg;
	else if (cmpstr(candidate, "lnorm")) *e = lnorm;
	else if (cmpstr(candidate, "jacobi")) *e = jacobi;
	else if (cmpstr(candidate, "spk")) *e = spk;
	else return 0;
	return 1;
}
int input(int argc, char*argv[], int *k, enum goal* e, double *** data, int *data_len, int *features_number) {
	FILE *f;

	if (argc != 4) return 0;

	*k = atoi(argv[1]);
	if (!getgoal(e, argv[2])) return 0;

	f = fopen(argv[3], "r");

	if (f == NULL) return 0;
	if ((*data_len = getdata(f, data, features_number)) == -1) return 0;
	fclose(f);
	return 1;
}

double **extractCentroids(double **data, int k) {
	int i, j;
	double **centroids;
	centroids = (double**)malloc(k * sizeof(double*));
	if (centroids == NULL) return NULL;
	for (i = 0; i < k; ++i) {
		centroids[i] = (double*)malloc(k * sizeof(double));
		if (centroids[i] == NULL)return NULL;
		for (j = 0; j < k; ++j) centroids[i][j] = data[i][j];
	}
	return centroids;
}

double ** calcFinalCentroids(double **data, int *map, int k, int data_length, int data_dim) {
	double ** sums;
	int i, j, *amounts;

	sums = (double**)malloc(k * sizeof(double*));
	amounts = (int*)calloc(k, sizeof(int));
	if (sums == NULL || amounts == NULL) return NULL;
	if (sums == NULL) return NULL;
	for (i = 0; i < k; ++i) {
		sums[i] = (double*)calloc(data_dim, sizeof(double));
		if (sums[i] == NULL)return NULL;
	}
	for (i = 0; i < data_length; ++i) {
		for (j = 0; j < data_dim; ++j) {
			sums[map[i]][j] += data[i][j];
		}
		amounts[map[i]] += 1;
	}
	for (i = 0; i < k; ++i) {
		for (j = 0; j < data_dim; ++j) {
			sums[i][j] /= amounts[i];
		}
	}
	return sums;
}

int main(int argc, char *argv[]) {
	int i, k, data_length, features, *original_to_Rnk_map;
	enum goal e;
	double *eigenvalues;
	double **eigenvectors, **initial_centroids, **final_centroids, 
		**data, **_wam, **_ddg, **_lnorm, **_Rnk;
	double  ***_jacobi;
	value_vector_map *vvmap;

	/*Set pointers to NULL, so free() method could
	 *ignore if a pointer was or was not used
	 */
	eigenvalues = NULL;
	eigenvectors = initial_centroids = final_centroids =
		data = _wam = _ddg = _lnorm = _Rnk = NULL;
	_jacobi = NULL;

	if (!input(argc, argv, &k, &e, &data, &data_length, &features)) {
		printf("Invalid Input!\n");
		return 0;
	}
	switch (e)
	{
	case wam:
		_wam = calcWAM(data, data_length, features);
		printArr(_wam, data_length, data_length);
		break;

	case ddg:
		_wam = calcWAM(data, data_length, features);
		_ddg = calcDDG(_wam, data_length);
		printArr(_ddg, data_length, data_length);

		break;
	case lnorm:
		_wam = calcWAM(data, data_length, features);
		_ddg = calcDDG(_wam, data_length);
		_lnorm = calcNGL(_wam, _ddg, data_length);
		printArr(_lnorm, data_length, data_length);

		break;
	case jacobi:
		/*Here data is of Mat(n,n) as an assumption*/
		_jacobi = calcJacobi(data, data_length);
		eigenvalues = _jacobi[0][0];
		eigenvectors = _jacobi[1];
		printArr(&eigenvalues, 1, data_length);
		printArr(eigenvectors, data_length, data_length);
		break;
	case spk:
		_wam = calcWAM(data, data_length, features);
		_ddg = calcDDG(_wam, data_length);
		_lnorm = calcNGL(_wam, _ddg, data_length);
		_jacobi = calcJacobi(_lnorm, data_length);
		eigenvalues = _jacobi[0][0];
		eigenvectors = _jacobi[1];
		vvmap = setMap(eigenvectors, eigenvalues, data_length);
		if (k == 0) k = determineK(vvmap, data_length);
		_Rnk = calcNormalaizedRnk(vvmap, data_length, k);
		initial_centroids = extractCentroids(_Rnk, k);
		original_to_Rnk_map = k_mean(_Rnk, initial_centroids, data_length, k, k, 300);
		final_centroids = calcFinalCentroids(data, original_to_Rnk_map, k, data_length, features);
		printArr(final_centroids, k, features);
		break;
	default:
		break;
	}
	
	
	/*free data*/
	for (i = 0; i < data_length; ++i)free(data[i]);
	free(data);
	/*free wam*/
	if (_wam != NULL) for (i = 0; i < data_length; ++i)free(_wam[i]);
	free(_wam);
	/*free ddg*/
	if (_ddg != NULL) for (i = 0; i < data_length; ++i)free(_ddg[i]);
	free(_ddg);
	/*free lnorm*/
	if (_lnorm != NULL) for (i = 0; i < data_length; ++i)free(_lnorm[i]);
	free(_ddg);
	/*free jacobi*/
	if (_jacobi != NULL) {
		free(_jacobi[0][0]); free(_jacobi[0]);
		for (i = 0; i < data_length; ++i)free(_jacobi[1][i]);
	}
	free(_jacobi);
	/*free vvmap*/
	free(vvmap);
	/*free Rnk*/
	if (_Rnk != NULL)for (i = 0; i < data_length; ++i)free(_Rnk[i]);
	free(_Rnk);
	/*free initial_centroids*/
	if (initial_centroids != NULL) for (i = 0; i < k; ++i) free(initial_centroids[i]);
	free(initial_centroids);
	/*free original_to_Rnk_map*/
	free(original_to_Rnk_map);
	/*free final_centroids*/
	if (final_centroids != NULL) for (i = 0; i < k; ++i) free(final_centroids[i]);
	free(final_centroids);
	return 0;
}