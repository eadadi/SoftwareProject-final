/*
This File Contains:
	1. Must-have includes:  for the c-api module
	2. Modularity includes: The required functions are implemented in external c files
	3. Python Objects to C Objects transformations: to make the usage of the c-api module available
	4. API-definitions:	defines the c-api behaviour
*/
//1.	Must-Have Includes
#define PY_SSIZE_T_CLEAN
#include<Python.h>
#include <stdbool.h>

//2.	Modularity Includes
#include "spkmeans.h"
#include "res/ddg.h"
#include "res/jacobi.h"
#include "res/lnorm.h"
#include "res/spk.h"
#include "res/wam.h"

/*		Remark.	Local Functions Declarations
static PyObject* c_ddg(PyObject *, PyObject *);
static PyObject* c_jacobi(PyObject *, PyObject *);
static PyObject* c_lnorm(PyObject *, PyObject *);
static PyObject* c_spk(PyObject *, PyObject *);
static PyObject* c_wam(PyObject *, PyObject *);*/

//3.0 Parse the Python Vectors to C Vectors:
void ** parsePythonArgs(PyObject *args) {
	int n, vectorDimension;
	double** c_vectors;
	void **arr;
	PyObject *vectors; Py_ssize_t i, j;
	
	arr = (void*)calloc(3, sizeof(void *));
	PyArg_ParseTuple(args, "Oii", &vectors, &n, &vectorDimension);

	c_vectors = (double**)calloc(n, sizeof(double*));
	for (i = 0; i < vectorDimension; ++i)
		c_vectors[i] = (double*)calloc(vectorDimension, sizeof(double));

	for (i = 0; i < n; ++i) { 	//build c_vectors arg
		PyObject *u = PyList_GetItem(vectors, i);
		for (j = 0; j < vectorDimension; ++j) {
			PyObject *value = PyList_GetItem(u, j);
			PyArg_Parse(value, "d", &c_vectors[i][j]);
		}
	}
	arr[0] = &c_vectors; arr[1] = &n, arr[2] = &vectorDimension;
	return arr;
}
void putPythonParsedArgsIntoVars(PyObject *args, double *** v, int *n, int*dim) {
	int *ptr1, *ptr2;
	double ***ptr0;
	void **parsedArr;
	parsedArr = parsePythonArgs(args);
	ptr0 = parsedArr[0]; ptr1 = parsedArr[1]; ptr2 = parsedArr[2];
	*v = *(ptr0); *n = *(ptr1); *dim = *(ptr2);
	free(parsedArr);
}

//3.1	Python Objects to C Objects Transformations
static PyObject* c_wam(PyObject *self, PyObject *args){
	double *row, **WAM;
	PyObject *pyWAM, *pyWAMrow;
	//Parse Python Args
	int n, vectorDimension, k, l;
	double **c_vectors;
	Py_ssize_t i, j;
	putPythonParsedArgsIntoVars(args, &c_vectors, &n, &vectorDimension);

	//calculate Weighted Adjacency Matrix
	WAM = calcWAM(c_vectors, n, vectorDimension);

	//define returned WAM Python Object
	pyWAM = PyList_New(n);
	for (i = 0, k = 0; i < n; ++i, ++k) {
		row = WAM[k];
		pyWAMrow = PyList_New(n);
		for (j = 0, l = 0; j < n; ++j, ++l) {
			PyObject *value;
			value = Py_BuildValue("d", row[l]);
			PyList_SetItem(pyWAMrow, j, value);
		}
		PyList_SetItem(pyWAM, i, pyWAMrow);
	}
	//free c_vectors, WAM
	for (k=0; i<n; ++k) {free(c_vectors[k]); free(WAM[k]);}
	free(c_vectors);free(WAM);
	return pyWAM;
}

static PyObject* c_ddg(PyObject *self, PyObject *args) {
	double *row, **WAM, **DDG;
	PyObject *pyDDG, *DDGrow;
	//Parse Python Args
	int n, vectorDimension, k, l;
	double **c_vectors;
	Py_ssize_t i, j;
	putPythonParsedArgsIntoVars(args, &c_vectors, &n, &vectorDimension);
	//calculate Weighted Adjacency Matrix
	WAM = calcWAM(c_vectors, n, vectorDimension);
	//calculate Diagonal Degree Matrix
	DDG = calcDDG(WAM, n);
	//define returned DDG Python Object
	pyDDG = PyList_New(n);
	for (i = 0, k = 0; i < n; ++i, ++k) {
		row = DDG[k];
		DDGrow = PyList_New(n);
		for (j = 0, l = 0; j < n; ++j, ++l) {
			PyObject *value;
			value = Py_BuildValue("d", row[l]);
			PyList_SetItem(DDGrow, j, value);
		}
		PyList_SetItem(pyDDG, i, DDGrow);
	}
	//free c_vectors, WAM, DDG
	for (k = 0; k < n; ++k) { free(c_vectors[k]); free(WAM[k]); free(DDG[k]); }
	free(c_vectors); free(WAM); free(DDG);
	return pyDDG;
}

static PyObject* c_lnorm(PyObject *self, PyObject *args){
	//Parse Python Args
	int n, vectorDimension, k, l;
	double **c_vectors;
	Py_ssize_t i, j;
	putPythonParsedArgsIntoVars(args, &c_vectors, &n, &vectorDimension);

	//calculate Weighted Adjacency Matrix
	double **WAM = calcWAM(c_vectors, n, vectorDimension);
	//calculate Diagonal Degree Matrix
	double **DDG = calcDDG(WAM, n);
	//calculate Normalized Graph Laplacian
	double **NGL = calcNGL(WAM, DDG, n);

	//define returned DDG Python Object
	PyObject* pyNGL = PyList_New(n);
	for (i = 0, k = 0; i < n; ++i, ++k) {
		double *row = NGL[k];
		PyObject *py_row = PyList_New(n);
		for (j = 0, l = 0; j < n; ++j, ++l) {
			PyObject *value;
			value = Py_BuildValue("d", row[l]);
			PyList_SetItem(py_row, j, value);
		}
		PyList_SetItem(pyNGL, i, py_row);
	}

	//free c_vectors, WAM, DDG, NGL
	for (k = 0; k < n; ++k) 
	{ free(c_vectors[k]); free(WAM[k]); free(DDG[k]); free(NGL[k]); }
	free(c_vectors); free(WAM); free(DDG); free(NGL);
	return pyNGL;
}

static PyObject* c_jacobi(PyObject *self, PyObject *args) {
	//Parse Python Args
	int n, vectorDimension, k, l;
	double **c_vectors;
	Py_ssize_t i, j;
	putPythonParsedArgsIntoVars(args, &c_vectors, &n, &vectorDimension);

	//calculate Weighted Adjacency Matrix
	double **WAM = calcWAM(c_vectors, n, vectorDimension);
	//calculate Diagonal Degree Matrix
	double **DDG = calcDDG(WAM, n);
	//calculate Normalized Graph Laplacian
	double **NGL = calcNGL(WAM, DDG, n);
	//do Jacobian Proccess
	double ***JACOBI = calcJacobi(NGL, n);
	double *eigenvalues = *JACOBI[0];
	double **eigenvectors = JACOBI[1];

	//define returned (eigenvalues,eigenvectors) Python Object
	PyObject* pyEIGENS = PyDict_New();
	//build python eigenvalues
	PyObject *py_eigenvalues = PyList_New(n);
	for (i = 0, k = 0; i < n; ++i, ++k) {
		PyObject *value;
		value = Py_BuildValue("d", eigenvalues[k]);
		PyList_SetItem(py_eigenvalues, i, value);
	}
	PyObject *key1 = Py_BuildValue("s", "eigenvalues");
	PyDict_SetItem(pyEIGENS, key1, py_eigenvalues);

	//build python eigenvectors
	PyObject *py_eigenvectors = PyList_New(n);
	for (i = 0, k = 0; i < n; ++i, ++k) {
		double *row = eigenvectors[k];
		PyObject *py_row = PyList_New(n);
		for (j = 0, l = 0; j < n; ++j, ++l) {
			PyObject *value;
			value = Py_BuildValue("d", row[l]);
			PyList_SetItem(py_row, j, value);
		}
		PyList_SetItem(py_eigenvectors, i, py_row);
	}
	PyObject *key2 = Py_BuildValue("s", "eigenvectors");
	PyDict_SetItem(pyEIGENS, key2, py_eigenvectors);

	//free c_vectors, WAM,, DDG, NGL
	for (k = 0; k < n; ++k)
	{
		free(c_vectors[k]); free(WAM[k]); free(DDG[k]); free(NGL[k]); 
		free(JACOBI[1][k]);
	}
	free(c_vectors); free(WAM); free(DDG); free(NGL); free(JACOBI[0][0]);
	free(JACOBI[0]); free(JACOBI[1]); free(JACOBI);
	return pyEIGENS;
}

static PyObject* c_spk(PyObject *self, PyObject *args){
	return Py_BuildValue("i",-1);
}

//4.	API-definitions
static PyMethodDef spkmeansMethods[] = {
	{"ddg",
	(PyCFunction)c_ddg,
	METH_VARARGS,
	PyDoc_STR("Diagonal Degree Matrix")},
	{"jacobi",
	(PyCFunction)c_jacobi,
	METH_VARARGS,
	PyDoc_STR("Eigenvalues and Eigenvectors")},
	{"lnorm",
	(PyCFunction)c_lnorm,
	METH_VARARGS,
	PyDoc_STR("Normalized Graph Lapliacan")},
	{"spk",
	(PyCFunction)c_spk,
	METH_VARARGS,
	PyDoc_STR("Full spectral kmeans")},
	{"wam",
	(PyCFunction)c_wam,
	METH_VARARGS,
	PyDoc_STR("Weighted Adjacency Matrix")},
	{NULL,NULL,0,NULL}
};

static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"spkmeans",
	NULL,
	-1,
	spkmeansMethods
};

PyMODINIT_FUNC
PyInit_spkmeans(void) {
	PyObject *m;
	m = PyModule_Create(&moduledef);
	if (!m) return NULL;
	return m;
}
