#define PY_SSIZE_T_CLEAN
#include<Python.h>
#include <stdbool.h>
/**
 * Python Stuff
 **/
static void k_mean(int , int ,double** , double ** , int , int );

static PyObject* fit(PyObject *self, PyObject *args){ 
	//Args from python:
	int k,max_iter,data_len,vector_len;
	PyObject *initial_centroids, *data;
	PyArg_ParseTuple(args, "iiiiOO", &k,&max_iter,&data_len,&vector_len,&initial_centroids,&data);
	//printf("k=%d max_iter=%d datalen=%d, vectorlen=%d\n",k,max_iter,data_len,vector_len);
	
	Py_ssize_t i, j;
	
	//Define C initial_centroids vector:
	double** c_initials = (double**)calloc(k,sizeof(double*));
	for (i=0; i<k; ++i){
		c_initials[i]= (double*)calloc(vector_len,sizeof(double));
	}
	//Build C initials_centroids vector
	for (i=0; i<k; ++i){
		PyObject *vector = PyList_GetItem(initial_centroids, i);
		if (!PyList_Check(vector)) printf("$$$ bug1\n");
		for (j=0; j<vector_len; ++j){
			PyObject* coordinate = PyList_GetItem(vector, j);
			PyArg_Parse(coordinate, "d", &c_initials[i][j]);
		}
	}

	//Define C data Matrice:
	double** c_data = (double**)calloc(data_len,sizeof(double*));
	for (i=0; i<data_len; ++i){
		c_data[i]= (double*)calloc(vector_len,sizeof(double));
	}
	//Build C data matrice:
	for (i=0; i<data_len; ++i){
		PyObject *vector = PyList_GetItem(data, i);
		if (!PyList_Check(vector)) printf("$$$ bug1\n");
		for (j=0; j<vector_len; ++j){
			PyObject* coordinate = PyList_GetItem(vector, j);
			PyArg_Parse(coordinate, "d", &c_data[i][j]);
		}
	}

	k_mean(k, max_iter, c_initials, c_data, data_len, vector_len);
	
	//Modify "initials" (which is pyobject) to hold the result:
	PyObject* result = PyList_New(PyList_Size(initial_centroids));
	if (!PyList_Check(result)) printf("bug2!!!\n");
	for (i=0; i<k; ++i){
		double* line = c_initials[i];
		PyObject *vector = PyList_New(vector_len);
		PyList_SetItem(result, i, vector);
		for (j=0; j<vector_len; ++j){
			PyObject* value;
			value = Py_BuildValue("d", line[j]);
			PyList_SetItem(vector, j, value);
		}
	}

	//free c_data & c_initials
	for (i=0; i<k; ++i) free(c_initials[i]);
	for (i=0; i<data_len; ++i) free(c_data[i]);
	free(c_initials);free(c_data);
    return result;
}

static PyMethodDef kmeansspMethods[] = {
	{"fit",
	(PyCFunction) fit,
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

/**
 * Code itself
 **/

static double distance(double * u, double * v, int dim) {
	double res; int i;
	res = 0;
	for (i=0;i<dim ; ++i) res+=(u[i]-v[i])*(u[i]-v[i]);
	return res;
}

static void add (double *** si, double *u,int * bound, int *cnt){
	int lim, i;
	double **b;
	if (*cnt == *bound) {
		lim = 2*(*bound);
		b = (double**)(calloc(lim,sizeof(double*)));
		assert(b!=NULL);
		for (i=0;i<*bound;++i){
			b[i] = (*si)[i];
		}
		free(*si);
		*si = b;
		*bound = lim;
	}
	(*si)[*cnt] = u;
	*cnt = *cnt+1;
}

static bool v_eq (double * u, double * v, int dim) {
	while(dim--) if (u[dim]!=v[dim]) return false;
	return true;
}

static double * calc_cent (double **si, int dim, int cnt, double * v){
	double * mean; int i, j;
	if(cnt>0){
		mean = (double*)calloc(dim,sizeof(double));
		assert(mean!=NULL);
		for (i=0; i<cnt; ++i){
			for (j=0;j<dim;++j) mean[j]+=si[i][j];
		}
		for (i=0;i<dim;++i) mean[i] = mean[i]/cnt;
		return mean;
	}
	else return v;
}

static void k_mean(int k, int max_iter,double** initials, double ** data, int mat_len, int dim){
	double** x = data;
	double **m = initials;
	double ***s, min_dis, temp, *temp_cent;
	int init, *bounds, *counts, argmin, i, j;
	bool changed =true;

	init = 128;
	changed = true;

	bounds = (int *)calloc(k,sizeof(int));
	assert(bounds!=NULL);
	counts = (int *)calloc(k,sizeof(int));
	assert(counts!=NULL);

	while (changed && max_iter ){
		changed = false;

		s = (double***)calloc(k,sizeof(double **));
		assert(s!=NULL);
		for (i = 0; i< k; i++) {
			s[i] = (double**)calloc(init,sizeof(double*));
			bounds[i] = init;
			counts[i]=0;
		}

		for (i=0;i<mat_len;++i){
			argmin = 0;
			min_dis = distance(x[i], m[0], dim);
			for (j=1 ; j<k; ++j){
				temp = distance(x[i], m[j] ,dim);
				if (temp < min_dis) {
					min_dis = temp;
					argmin = j;
				}
			}
			add(&s[argmin],x[i], &bounds[argmin], &counts[argmin]);
		}
		for (i=0; i<k; ++i){
			temp_cent = calc_cent(s[i], dim, counts[i],m[i]);
			if (!v_eq(temp_cent,m[i],dim)) {
				changed = true;
				for(j=0;j<dim;++j) m[i][j] = temp_cent[j];
			}
			free(temp_cent);
		}
		max_iter-=1;

		for (i=0;i<k; ++i) free(s[i]);
		free(s);
	}
	
	

	free(counts);free(bounds);
	//for(i=0;i<k;++i) free(m[i]);
	//free(m);
	//for (i=0;i<mat_len; ++i) free(x[i]);
	
	//Result(m, k, dim);
	//return m;
}

