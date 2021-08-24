#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

static double ** static_py_calcWAM(double ** vectors, int n, int vectorDimension) {
	return calcWAM(vectors, n, vectorDimension);
}
static double **static_py_calcDDG(double ** M, int n) {
	return calcDDG(M, n);
}
static double ** static_py_calcNGL(double ** W, double ** D, int n) {
	return calcNGL(W, D, n);
}
static double*** static_py_calcJacobi(double **A, int n) {
	return calcJacobi(A, n);
}
static value_vector_map* static_py_setMap(double ** eigenvectors, double * eigenvalues, int n) {
	return setMap(eigenvectors, eigenvalues, n);
}
static int static_py_determineK(value_vector_map *map, int n) {
	return determineK(map, n);
}
static double ** static_py_calcNormalaizedRnk(value_vector_map *map, int n, int clustersNumber) {
	return calcNormalaizedRnk(map, n, clustersNumber);
}
int* static_py_k_mean(double **datapoints, double **centroids, int datapoints_amount,
	int clusters_amount, int datapoint_length, int max_iter) {
	return  k_mean(datapoints, centroids, datapoints_amount,	
				clusters_amount, datapoint_length, max_iter);
}

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
	*datapoints_amount = (int) _datapoints_amount;
	_clustersAmount = PyList_Size(_Centroids);
	*clusters_amount = (int) _clustersAmount;
	_vector = PyList_GetItem(_Datapoints, 0);
	_vectorDim = PyList_Size(_vector);
	*datapoint_length = (int) _vectorDim;


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
	WAM = static_py_calcWAM(vectors, n, vectorDimension);
	
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
	WAM = static_py_calcWAM(vectors, n, vectorDimension);
	if (WAM == NULL) return NULL;
	DDG = static_py_calcDDG(WAM, n);
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
	WAM = static_py_calcWAM(vectors, n, vectorDimension);
	if (WAM == NULL) return NULL;
	DDG = static_py_calcDDG(WAM, n);
	if (DDG == NULL) return NULL;
	NGL = static_py_calcNGL(WAM, DDG, n);
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
	JACOBI = static_py_calcJacobi(vectors,n);
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
	PyObject *_number, *_RNK, *_RNKrow, *_vvmap;
	Py_ssize_t i, j, _n, _clusterNumber;
	double **vectors, **WAM, **DDG, **NGL,
		***JACOBI, *eigenvalues, **eigenvectors, **RNK;
	value_vector_map *vvmap;
	vectors = NULL;

	if (EasyArgParse(args, &vectors, &n, &vectorDimension, &clustersNumber) == 0) return NULL;
	_n = PyLong_AsSize_t(Py_BuildValue("l", n));
	
	// Computation
	WAM = static_py_calcWAM(vectors, n, vectorDimension);
	if (WAM == NULL) return NULL;
	//Freeing vectors
	for (k = 0; k < n; ++k) free(vectors[k]);
	free(vectors);

	DDG = static_py_calcDDG(WAM, n);
	if (DDG == NULL) return NULL;
	NGL = static_py_calcNGL(WAM, DDG, n);
	if (NGL == NULL) return NULL;
	//Freeing DDG
	for (k = 0; k < n; ++k) free(DDG[k]);
	free(DDG);
	//Freeing WAM
	for (k = 0; k < n; ++k) free(WAM[k]);
	free(WAM);

	JACOBI = static_py_calcJacobi(NGL, n);
	if (JACOBI == NULL) return NULL;
	eigenvalues = JACOBI[0][0];
	eigenvectors = JACOBI[1];

	vvmap = static_py_setMap(eigenvectors, eigenvalues, n);
	if (vvmap == NULL) return NULL;

	if (clustersNumber == 0) clustersNumber = static_py_determineK(vvmap, n);
	_clusterNumber = PyLong_AsSize_t(Py_BuildValue("l", clustersNumber));

	RNK = static_py_calcNormalaizedRnk(vvmap, n, clustersNumber);
	if (RNK == NULL) return NULL;
	
	//Parsing RNK to Python
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

	//Parsing vvmap (indexes only) to Python
	_vvmap = PyList_New(_n);
	if (!PyList_Check(_RNK)) return NULL;
	for (i = 0; i < _n; ++i) {
		_number = Py_BuildValue("i", vvmap[i].index);
		if (!_number) return NULL;
		if (PyList_SetItem(_vvmap, i, _number) == -1) return NULL;
	}
	//Freeing NormalaizedRnk
	for (k = 0; k < n; ++k) free(RNK[k]);
	free(RNK);
	//Freeing value_vector_map
	free(vvmap);
	//Freeing JACOBI[0]
	free(JACOBI[0][0]); free(JACOBI[0]);
	//Freeing JACOBI[1]
	for (k = 0; k < n; ++k) free(JACOBI[1][k]); free(JACOBI[1]);
	free(JACOBI);
	//Freeing NGL
	for (k = 0; k < n; ++k) free(NGL[k]);
	free(NGL);
	
	//Create (RNK,vvmap) Tuple:
	//return _RNK;
	return Py_BuildValue("(OO)", _RNK, _vvmap);
}
static PyObject* c_kmeans(PyObject *self, PyObject *args) {
	int datapoints_amount, clusters_amount, datapoint_length, max_iter, succ, *map;
	PyObject *_MAP, *_number;
	double **datapoints, **centroids;
	succ = kmeans_EasyArgParse(args, &datapoints, &centroids, &datapoints_amount,
		&clusters_amount, &datapoint_length, &max_iter);
	Py_ssize_t i, _datapointsAmount;
	if (succ == 0) return NULL;
	
	_datapointsAmount = PyLong_AsSize_t(Py_BuildValue("n",datapoints_amount));

	map = static_py_k_mean(datapoints, centroids, datapoints_amount,
		clusters_amount, datapoint_length, max_iter);

	//Parsing to Python
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


