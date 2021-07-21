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

//3.	Python Objects to C Objects Transformations

static PyObject* c_wam(PyObject *self, PyObject *args){
	int n, vectorDimension;
	PyObject *vectors;
	double** c_vectors;
	Py_ssize_t i,j;

	PyArg_ParseTuple(args, "Oii", &vectors, &n, &vectorDimension);
	
	//define c_vectors arg
	c_vectors = (double**)calloc(n,sizeof(double*));
	for (i=0; i<vectorDimension; ++i) 
	c_vectors[i] = (double*)calloc(vectorDimension,sizeof(double));

	//build c_vectors arg
	
	for (i=0; i<n; ++i) {
		PyObject *u = PyList_GetItem(vectors, i);
		for (j=0;j<vectorDimension;++j) {
			PyObject *value = PyList_GetItem(u, j);
			PyArg_Parse(value, "d", &c_vectors[i][j]);
		}
	}

	//calculate Weighted Adjacency Matrix
	double **WAM = calcWAM(c_vectors, n, vectorDimension);

	//define returnedWAM Python Object
	PyObject* pyWAM = PyList_New(n);
	for (i=0; i<n;++i){
		double *row = WAM[i];
		PyObject *py_wam_row = PyList_New(n);
		for (j=0;j<n;++j){
			PyObject *value;
			value = Py_BuildValue("d", row[j]);
			PyList_SetItem(py_wam_row, j, value);
		}
		PyList_SetItem(pyWAM, i, py_wam_row);
	}

	//free c_vectors, WAM
	for (int i=0; i<n; ++i) {free(c_vectors[i]); free(WAM[i]);}
	free(c_vectors);free(WAM);
	return pyWAM;
}

static PyObject* c_ddg(PyObject *self, PyObject *args) {
	return Py_BuildValue("i",-1);
}

static PyObject* c_jacobi(PyObject *self, PyObject *args){
	return Py_BuildValue("i",-1);
}

static PyObject* c_lnorm(PyObject *self, PyObject *args){
	return Py_BuildValue("i",-1);
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
