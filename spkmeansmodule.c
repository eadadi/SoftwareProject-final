#define PY_SSIZE_T_CLEAN
#include <Python.h>


#include"res/tools.c"
#include "res/value_vector_map_struct.h"
#include "spkmeans.h"
#include "res/wam.c"
#include "res/ddg.c"
#include "res/lnorm.c"
#include "res/jacobi.c"
#include "res/eigenpap.c"
#include "res/spk.c"
#include "res/fit.c"
/*
static double distance(double *, double *, int);
static double **calcWAM(double **, int, int);
static double **calcDDG(double ** M, int n);
static double ** raise_D_diagonal_in_minus_half(double **D, int n);
static double ** left_multip_of_diagonal_matrice(double ** D, double **M, int n);
static double ** right_multip_of_diagonal_matrice(double ** M, double **D, int n);
static double ** calcNGL(double ** W, double ** D, int n);
#define EPSILON 0.001
#define JACOBI_MAX_ITERATIONS_NUMBER 100
static double*** calcJacobi(double **A, int n);
static double ** buildPij(double **M, int n, int i, int j, double c, double s);
static int* getPivotIndexes(double ** M, int n);
static double* calcCandS(double **M, int i, int j);
static double ** MatriceMultiplication(double ** A, double **B, int n);
static double ** calcAtag(double ** A, int n, int i, int j, double c, double s);
static int isDiagonalEnough(double **A, double **B, int n);
static int cmpfunc(const void *a, const void *b);
static void sortMap(value_vector_map *map, int n);
static int determineK(value_vector_map *map, int n);
static double ** calcNormalaizedRnk(value_vector_map *map, int n, int clustersNumber);
static int* k_mean(double **datapoints, double **centroids, int datapoints_amount,
	int clusters_amount, int datapoint_length, int max_iter);
static int getClosestCluster(double *candidate, double**centroids, int clusters_amount,
	int datapoint_length);
static int updateCentroids(double **datapoints, double **centroids,
	int datapoints_amount, int clusters_amount, int *data_to_centroids_map, int datapoint_length);
static void add_u_to_v(double *u, double *v, int dim);
*/


static int EasyArgParse(PyObject *args, double ***Vectors, int *num_of_vectors, int *vectorDim, int *k_argument) {
	int n, vectorDimension, k, l;
	PyObject *_vectors, *_vector, *_number;
	Py_ssize_t i, j, _n, _vectorDim;
	double **vectors, *vector;

	if (!PyArg_ParseTuple(args, "Oiii", &_vectors, &n, &vectorDimension, &k)) return 0;
	if (!PyList_Check(_vectors)) return 0;
	_n = PyList_Size(_vectors);

	vectors = (double**)calloc(n, sizeof(double*));
	if (vectors == NULL) return 0;
	for (l = 0; l < n; ++l) {
		vectors[l] = (double*)calloc(vectorDimension, sizeof(double));
		if(vectors[l] == NULL) return 0;
	}

	for (i = 0; i < _n; i++) {
		_vector = PyList_GetItem(_vectors, i);
		if (!PyList_Check(_vector)) return 0;
		_vectorDim = PyList_Size(_vector);
		vector = vectors[i];
		for (j = 0; j < _vectorDim; ++j) {
			_number = PyList_GetItem(_vector, j);
			vector[j] = PyFloat_AsDouble(_number);
			if (vector[j] == -1 && PyErr_Occurred()) return 0;
		}
	}
	*Vectors = vectors; *num_of_vectors = n; *vectorDim = vectorDimension;
	if (k_argument != NULL) *k_argument = k;
	return 1;
}
static int kmeans_EasyArgParse(PyObject *args, double ***Datapoints, double ***Centroids,
	int *datapoints_amount, int *clusters_amount, int *datapoint_length, int *max_iter) {
	int l;
	PyObject *_Datapoints, *_Centroids, *_vector, *_number;
	Py_ssize_t i, j, _datapoints_amount, _vectorDim, _clustersAmount;
	double **vectors, **centros, *vector;
	if (!PyArg_ParseTuple(args, "OOi", &_Datapoints,&_Centroids, max_iter)) return 0;
	if (!PyList_Check(_Datapoints)) return 0;
	if (!PyList_Check(_Centroids)) return 0;
	_datapoints_amount = PyList_Size(_Datapoints);
	*datapoints_amount = _datapoints_amount;
	_clustersAmount = PyList_Size(_Centroids);
	*clusters_amount = _clustersAmount;
	_vector = PyList_GetItem(_Datapoints, 0);
	_vectorDim = PyList_Size(_vector);
	*datapoint_length = _vectorDim;


	/*
	Datapoints Parsing
	*/
	vectors = (double**)calloc(*datapoints_amount, sizeof(double*));
	if (vectors == NULL) return 0;
	for (l = 0; l < *datapoints_amount; ++l) {
		vectors[l] = (double*)calloc(*datapoint_length, sizeof(double));
		if (vectors[l] == NULL) return 0;
	}

	for (i = 0; i < _datapoints_amount; i++) {
		_vector = PyList_GetItem(_Datapoints, i);
		if (!PyList_Check(_vector)) return 0;
		_vectorDim = PyList_Size(_vector);
		vector = vectors[i];
		for (j = 0; j < _vectorDim; ++j) {
			_number = PyList_GetItem(_vector, j);
			vector[j] = PyFloat_AsDouble(_number);
			if (vector[j] == -1 && PyErr_Occurred()) return 0;
		}
	}
	/*
	Centroids Parsing
	*/
	centros = (double**)calloc(*clusters_amount, sizeof(double*));
	if (centros == NULL) return 0;
	for (l = 0; l < *clusters_amount; ++l) {
		centros[l] = (double*)calloc(*datapoint_length, sizeof(double));
		if (centros[l] == NULL) return 0;
	}

	for (i = 0; i < _clustersAmount; i++) {
		_vector = PyList_GetItem(_Centroids, i);
		if (!PyList_Check(_vector)) return 0;
		_vectorDim = PyList_Size(_vector);
		vector = centros[i];
		for (j = 0; j < _vectorDim; ++j) {
			_number = PyList_GetItem(_vector, j);
			vector[j] = PyFloat_AsDouble(_number);
			if (vector[j] == -1 && PyErr_Occurred()) return 0;
		}
	}

	*Datapoints = vectors; *Centroids = centros;
	return 1;
}
static PyObject* c_wam(PyObject *self, PyObject *args) {
	int n, vectorDimension, k;
	PyObject *_number, *_WAM, *_WAMrow;
	Py_ssize_t i, j, _n;
	double **vectors, **WAM;
	vectors = NULL;

	if (EasyArgParse(args, &vectors, &n, &vectorDimension, NULL) == 0) return NULL;
	_n = PyLong_AsSize_t(Py_BuildValue("l", n));
	// Computation
	WAM = calcWAM(vectors, n, vectorDimension);
	
	//Parsing to Python
	_WAM = PyList_New(_n);
	if (!PyList_Check(_WAM)) return NULL;
	for (i = 0; i < _n; ++i) {
		_WAMrow = PyList_New(_n);
		if (!PyList_Check(_WAMrow)) return NULL;
		for (j = 0; j < _n; ++j) {
			_number = Py_BuildValue("d", WAM[i][j]);
			if (!_number) return NULL;
			PyList_SetItem(_WAMrow, j, _number);
		}
		PyList_SetItem(_WAM, i, _WAMrow);
	}

	//Freeing WAM
	for (k = 0; k < n; ++k) free(WAM[k]);
	free(WAM);
	//Freeing vectors
	for (k = 0; k < n; ++k) free(vectors[k]);
	free(vectors);
	return _WAM;
}
static PyObject* c_ddg(PyObject *self, PyObject *args) {
	int n, vectorDimension, k;
	PyObject *_number, *_DDG, *_DDGrow;
	Py_ssize_t i, j, _n;
	double **vectors, **WAM, **DDG;
	vectors = NULL;

	if (EasyArgParse(args, &vectors, &n, &vectorDimension, NULL) == 0) return NULL;
	_n = PyLong_AsSize_t(Py_BuildValue("l", n));
	// Computation
	WAM = calcWAM(vectors, n, vectorDimension);
	if (WAM == NULL) return NULL;
	DDG = calcDDG(WAM, n);
	if (DDG == NULL) return NULL;
	//Parsing to Python
	_DDG = PyList_New(_n);
	if (!PyList_Check(_DDG)) return NULL;
	for (i = 0; i < _n; ++i) {
		_DDGrow = PyList_New(_n);
		if (!PyList_Check(_DDGrow)) return NULL;
		for (j = 0; j < _n; ++j) {
			_number = Py_BuildValue("d", DDG[i][j]);
			if (!_number) return NULL;
			if(PyList_SetItem(_DDGrow, j, _number) ==-1) return NULL;
		}
		if(PyList_SetItem(_DDG, i, _DDGrow)==-1) return NULL;
	}
	//Freeing DDG
	for (k = 0; k < n; ++k) free(DDG[k]);
	free(DDG);
	//Freeing WAM
	for (k = 0; k < n; ++k) free(WAM[k]);
	free(WAM);
	//Freeing vectors
	for (k = 0; k < n; ++k) free(vectors[k]);
	free(vectors);
	return _DDG;
}
static PyObject* c_lnorm(PyObject *self, PyObject *args) {
	int n, vectorDimension, k;
	PyObject *_number, *_NGL, *_NGLrow;
	Py_ssize_t i, j, _n;
	double **vectors, **WAM, **DDG, **NGL;
	vectors = NULL;

	if (EasyArgParse(args, &vectors, &n, &vectorDimension, NULL) == 0) return NULL;
	_n = PyLong_AsSize_t(Py_BuildValue("l", n));
	// Computation
	WAM = calcWAM(vectors, n, vectorDimension);
	if (WAM == NULL) return NULL;
	DDG = calcDDG(WAM, n);
	if (DDG == NULL) return NULL;
	NGL = calcNGL(WAM, DDG, n);
	if (NGL == NULL) return NULL;

	//Parsing to Python
	_NGL = PyList_New(_n);
	if (!PyList_Check(_NGL)) return NULL;
	for (i = 0; i < _n; ++i) {
		_NGLrow = PyList_New(_n);
		if (!PyList_Check(_NGLrow)) return NULL;
		for (j = 0; j < _n; ++j) {
			_number = Py_BuildValue("d", NGL[i][j]);
			if (!_number) return NULL;
			if (PyList_SetItem(_NGLrow, j, _number) == -1) return NULL;
		}
		if (PyList_SetItem(_NGL, i, _NGLrow) == -1) return NULL;
	}
	//Freeing NGL
	for (k = 0; k < n; ++k) free(NGL[k]);
	free(NGL);
	//Freeing DDG
	for (k = 0; k < n; ++k) free(DDG[k]);
	free(DDG);
	//Freeing WAM
	for (k = 0; k < n; ++k) free(WAM[k]);
	free(WAM);
	//Freeing vectors
	for (k = 0; k < n; ++k) free(vectors[k]);
	free(vectors);
	return _NGL;
}
static PyObject* c_jacobi(PyObject *self, PyObject *args) {
	int n, vectorDimension, k;
	PyObject *_number,*_JACOBI,
		*_eigenvalues, *_eigenvectors, *_eigenvectorsRow;
	Py_ssize_t i, j, _n;
	double **vectors,  ***JACOBI, *eigenvalues, **eigenvectors;
	vectors = NULL;

	if (EasyArgParse(args, &vectors, &n, &vectorDimension, NULL) == 0) return NULL;
	_n = PyLong_AsSize_t(Py_BuildValue("l", n));
	/**
	*Note that in this case vectors dimension is 'n' because we got n*n Matrice
	*/
	// Computation
	JACOBI = calcJacobi(vectors,n);
	if (JACOBI == NULL) return NULL;
	eigenvalues = JACOBI[0][0];
	eigenvectors = JACOBI[1];
	
	//Parsing eigenvectors to Python
	_eigenvectors = PyList_New(_n);
	if (!PyList_Check(_eigenvectors)) return NULL;
	for (i = 0; i < _n; ++i) {
		_eigenvectorsRow = PyList_New(_n);
		if (!PyList_Check(_eigenvectorsRow)) return NULL;
		for (j = 0; j < _n; ++j) {
			_number = Py_BuildValue("d", eigenvectors[i][j]);
			if (!_number) return NULL;
			if (PyList_SetItem(_eigenvectorsRow, j, _number) == -1) return NULL;
		}
		if (PyList_SetItem(_eigenvectors, i, _eigenvectorsRow) == -1) return NULL;
	}
	//Parsing eigenvalues to Python
	_eigenvalues = PyList_New(_n);
	if (!PyList_Check(_eigenvalues)) return NULL;
	for (j = 0; j < _n; ++j) {
		_number = Py_BuildValue("d", eigenvalues[j]);
		if (!_number) return NULL;
		if (PyList_SetItem(_eigenvalues, j, _number) == -1) return NULL;
	}
	//Combine together
	_JACOBI = PyList_New(2);
	if (!PyList_Check(_JACOBI)) return NULL;
	if (PyList_SetItem(_JACOBI, 0, _eigenvalues) == -1) return NULL;
	if (PyList_SetItem(_JACOBI, 1, _eigenvectors) == -1) return NULL;


	//Freeing JACOBI[0]
	free(JACOBI[0][0]); free(JACOBI[0]);
	//Freeing JACOBI[1]
	for (k = 0; k < n; ++k) free(JACOBI[1][k]); free(JACOBI[1]);
	free(JACOBI);
	//Freeing vectors
	for (k = 0; k < n; ++k) free(vectors[k]);
	free(vectors);
	return _JACOBI;
}
static PyObject* c_spk(PyObject *self, PyObject *args) {
	int n, vectorDimension, k, clustersNumber;
	PyObject *_number, *_RNK, *_RNKrow;
	Py_ssize_t i, j, _n, _clusterNumber;
	double **vectors, **WAM, **DDG, **NGL,
		***JACOBI, *eigenvalues, **eigenvectors, **RNK;
	value_vector_map *map;
	vectors = NULL;

	if (EasyArgParse(args, &vectors, &n, &vectorDimension, &clustersNumber) == 0) return NULL;
	_n = PyLong_AsSize_t(Py_BuildValue("l", n));
	
	// Computation
	WAM = calcWAM(vectors, n, vectorDimension);
	if (WAM == NULL) return NULL;
	//Freeing vectors
	for (k = 0; k < n; ++k) free(vectors[k]);
	free(vectors);

	DDG = calcDDG(WAM, n);
	if (DDG == NULL) return NULL;
	NGL = calcNGL(WAM, DDG, n);
	if (NGL == NULL) return NULL;
	//Freeing DDG
	for (k = 0; k < n; ++k) free(DDG[k]);
	free(DDG);
	//Freeing WAM
	for (k = 0; k < n; ++k) free(WAM[k]);
	free(WAM);

	JACOBI = calcJacobi(NGL, n);
	if (JACOBI == NULL) return NULL;
	eigenvalues = JACOBI[0][0];
	eigenvectors = JACOBI[1];

	map = setMap(eigenvectors, eigenvalues, n);
	if (map == NULL) return NULL;

	if (clustersNumber == 0) clustersNumber = determineK(map, n);
	_clusterNumber = PyLong_AsSize_t(Py_BuildValue("l", clustersNumber));

	RNK = calcNormalaizedRnk(map, n, clustersNumber);
	if (RNK == NULL) return NULL;
	
	//Parsing to Python
	_RNK = PyList_New(_n);
	if (!PyList_Check(_RNK)) return NULL;
	for (i = 0; i < _n; ++i) {
		_RNKrow = PyList_New(_clusterNumber);
		if (!PyList_Check(_RNKrow)) return NULL;
		for (j = 0; j < _clusterNumber; ++j) {
			_number = Py_BuildValue("d", RNK[i][j]);
			if (!_number) return NULL;
			if (PyList_SetItem(_RNKrow, j, _number) == -1) return NULL;
		}
		if (PyList_SetItem(_RNK, i, _RNKrow) == -1) return NULL;
	}


	//Freeing NormalaizedRnk
	for (k = 0; k < n; ++k) free(RNK[k]);
	free(RNK);
	//Freeing value_vector_map
	free(map);
	//Freeing JACOBI[0]
	free(JACOBI[0][0]); free(JACOBI[0]);
	//Freeing JACOBI[1]
	for (k = 0; k < n; ++k) free(JACOBI[1][k]); free(JACOBI[1]);
	free(JACOBI);
	//Freeing NGL
	for (k = 0; k < n; ++k) free(NGL[k]);
	free(NGL);
	
	return _RNK;
}
static PyObject* c_kmeans(PyObject *self, PyObject *args) {
	int succ;
	int datapoints_amount, clusters_amount, datapoint_length, max_iter;
	PyObject *_MAP, *_number;
	double **datapoints, **centroids;
	succ = kmeans_EasyArgParse(args, &datapoints, &centroids, &datapoints_amount,
		&clusters_amount, &datapoint_length, &max_iter);
	Py_ssize_t i, k, _clustersAmount, _datapointsAmount, _datapointLength;
	int *map;
	if (succ == 0) return NULL;
	
	_datapointsAmount = PyLong_AsSize_t(Py_BuildValue("n",datapoints_amount));
	_clustersAmount = PyLong_AsSize_t(Py_BuildValue("n", clusters_amount));
	_datapointLength = PyLong_AsSize_t(Py_BuildValue("n", datapoint_length));


	map = k_mean(datapoints, centroids, datapoints_amount, clusters_amount, datapoint_length, max_iter);

	//Parsing to Python
	k = 1;
	_MAP = PyList_New(_datapointsAmount);
	if (!PyList_Check(_MAP)) return NULL;
	for (i = 0; i < _datapointsAmount; ++i) {
		_number = Py_BuildValue("l", map[i]);
		if (!_number) return NULL;
		if (PyList_SetItem(_MAP, i, _number) == -1) return NULL;
	}

	//free map
	free(map);
	//free datapoints
	for (i = 0; i < datapoints_amount; ++i) free(datapoints[i]);
	free(datapoints);
	//free centroids
	for (i = 0; i < clusters_amount; ++i) free(centroids[i]);
	free(centroids);
	return _MAP;
}

static PyMethodDef spkmeansMethods[] = {
	{"ddg",(PyCFunction)c_ddg,METH_VARARGS,PyDoc_STR("Diagonal Degree Matrix")},
	{"jacobi",(PyCFunction)c_jacobi,METH_VARARGS,PyDoc_STR("Eigenvalues and Eigenvectors")},
	{"lnorm",(PyCFunction)c_lnorm,METH_VARARGS,PyDoc_STR("Normalized Graph Lapliacan")},
	{"spk",(PyCFunction)c_spk,METH_VARARGS,PyDoc_STR("Full spectral kmeans")},
	{"wam",(PyCFunction)c_wam,METH_VARARGS,PyDoc_STR("Weighted Adjacency Matrix")},
	{"kmeans",(PyCFunction)c_kmeans,METH_VARARGS,PyDoc_STR("Kmeans Alg")},
	{NULL,NULL,0,NULL}
};

static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,"spkmeans",NULL,-1,spkmeansMethods
};

PyMODINIT_FUNC
PyInit_spkmeans(void) {
	PyObject *m;
	m = PyModule_Create(&moduledef);
	if (!m) return NULL;
	return m;
}


