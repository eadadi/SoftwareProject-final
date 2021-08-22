#include "spkmeans.h"

#include "res/tools.c"
#include "res/wam.c"
#include "res/ddg.c"
#include "res/lnorm.c"
#include "res/jacobi.c"
#include "res/spk.c"
#include "res/eigenpap.c"
#include "res/fit.c"
#include "res/c_interface.c"

void printArr(double ** A, int n, int m) {
	int i, j;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < m; ++j) {
			printf("%.4f", A[i][j]);
			j + 1 == m ? printf("\n") : printf(",");
		}
	}
}

int main(int argc, char *argv[]) {
	int i, k, data_length, features, *original_to_Rnk_map, error;
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
	original_to_Rnk_map = NULL;

	if (!input(argc, argv, &k, &e, &data, &data_length, &features)) {
		printf("Invalid Input!\n");
		return 0;
	}
	error = 0;
	switch (e)
	{
	case wam:
		_wam = calcWAM(data, data_length, features); 
		if (_wam == NULL) { error = 1; break;}
		printArr(_wam, data_length, data_length);
		break;

	case ddg:
		_wam = calcWAM(data, data_length, features);
		if (_wam == NULL) { error = 1; break; }
		_ddg = calcDDG(_wam, data_length);
		if (_ddg == NULL) { error = 1; break; }
		printArr(_ddg, data_length, data_length);
		break;
	case lnorm:
		_wam = calcWAM(data, data_length, features);
		if (_wam == NULL) { error = 1; break; }
		_ddg = calcDDG(_wam, data_length);
		if (_ddg == NULL) { error = 1; break; }
		_lnorm = calcNGL(_wam, _ddg, data_length);
		if (_lnorm == NULL) { error = 1; break; }
		printArr(_lnorm, data_length, data_length);

		break;
	case jacobi:
		/*Here data is of Mat(n,n) as an assumption*/
		_jacobi = calcJacobi(data, data_length);
		if (_jacobi == NULL) { error = 1; break; }
		eigenvalues = _jacobi[0][0];
		eigenvectors = _jacobi[1];
		printArr(eigenvectors, data_length, data_length);
		printArr(&eigenvalues, 1, data_length);
		break;
	case spk:
		_wam = calcWAM(data, data_length, features);
		if (_wam == NULL) { error = 1; break; }
		_ddg = calcDDG(_wam, data_length);
		if (_ddg == NULL) { error = 1; break; }
		_lnorm = calcNGL(_wam, _ddg, data_length);
		if (_lnorm == NULL) { error = 1; break; }
		_jacobi = calcJacobi(_lnorm, data_length);
		if (_jacobi == NULL) { error = 1; break; }
		eigenvalues = _jacobi[0][0];
		eigenvectors = _jacobi[1];
		vvmap = setMap(eigenvectors, eigenvalues, data_length);
		if (vvmap == NULL) { error = 1; break; }
		if (k == 0) k = determineK(vvmap, data_length);
		_Rnk = calcNormalaizedRnk(vvmap, data_length, k);
		if (_Rnk == NULL) { error = 1; break; }
		initial_centroids = extractCentroids(_Rnk, k);
		if (initial_centroids == NULL) { error = 1; break; }
		original_to_Rnk_map = k_mean(_Rnk, initial_centroids, data_length, k, k, 300);
		if (original_to_Rnk_map == NULL) { error = 1; break; }
		final_centroids = calcFinalCentroids(data, original_to_Rnk_map, k, data_length, features);
		if (final_centroids == NULL) { error = 1; break; }
		printArr(final_centroids, k, features);
		break;
	default:
		error = 1;
		break;
	}
	
	if (error == 1) { printf("An Error Has Occured\n"); return 0; }

	/*free data*/
	if (data != NULL)for (i = 0; i < data_length; ++i)free(data[i]);
	free(data);
	/*free wam*/
	if (_wam != NULL) for (i = 0; i < data_length; ++i)free(_wam[i]);
	free(_wam);
	/*free ddg*/
	if (_ddg != NULL) for (i = 0; i < data_length; ++i)free(_ddg[i]);
	free(_ddg);
	/*free lnorm*/
	if (_lnorm != NULL) for (i = 0; i < data_length; ++i)free(_lnorm[i]);
	free(_lnorm);
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
			/*sums[map[i]][j] += data[i][j];*/
		}
		amounts[map[i]] += 1;
	}
	for (i = 0; i < k; ++i) {
		for (j = 0; j < data_dim; ++j) {
			if (amounts[i] > 0)sums[i][j] /= amounts[i];
		}
	}

	free(amounts);
	return sums;
}