/*
This File Contains:
	1. Must-have includes:  for the c-api module
	2. Modularity includes: The required functions are implemented in external c files
	3. Function Declarations
	4. API-definitions:	defines the c-api behaviour
	5. Python Objects to C Objects transformations: to make the usage of the c-api module available
*/
//1.	Must-Have Includes
#define PY_SSIZE_T_CLEAN
#include<Python.h>
#include <stdbool.h>

//2.	Modularity Includes


//3.	Functions Declarations


//4.	API-definitions
static PyMethodDef kmeansspMethods[] = {
	{"fit",
	(PyCFunction)fit,
	METH_VARARGS,
	PyDoc_STR("The link between my python code and c code EA15")},
	{NULL,NULL,0,NULL}
};

static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"mykmeanssp",
	NULL,
	-1,
	kmeansspMethods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void) {
	PyObject *m;
	m = PyModule_Create(&moduledef);
	if (!m) return NULL;
	return m;
}

//5.	Python Objects to C Objects Transformations