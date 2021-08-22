#include "../spkmeans.h"

int parent_loc(int i) {
	int l;
	l = i / 2;
	return l;
}


double parent(matrice_max_heap *h, int i) {
	int l;
	l = parent_loc(i);
	return h->values[l];
}


int left_loc(int i) {
	int l;
	l = 2 * i;
	return l;
}


double left(matrice_max_heap *h, int i) {
	int l;
	l = left_loc(i);
	return h->values[l];
}


int right_loc(int i) {
	int l;
	l = 2 * i + 1;
	return l;
}


double right(matrice_max_heap *h, int i) {
	int l;
	l = right_loc(i);
	return h->values[l];
}


double key(matrice_max_heap *h, int i) {
	return h->values[i];
}


void set_key(matrice_max_heap *h, int i, double k) {
	h->values[i] = k;
}


void set_indexes(matrice_max_heap *h, int index, int i, int j) {
	h->values_to_mat[index * HEAP_MEM] = i;
	h->values_to_mat[index * HEAP_MEM + 1] = j;
}


int heap_len(matrice_max_heap *h) {
	return h->values_to_mat[0];
}


void interswitch(matrice_max_heap *h, int i, int j) {
	int j_ind0, j_ind1, i_ind0, i_ind1, tmp;
	double j_k, i_k;
	j_k = key(h, j);
	j_ind0 = h->values_to_mat[HEAP_MEM * j];
	j_ind1 = h->values_to_mat[HEAP_MEM * j + 1];

	i_k = key(h, i);
	i_ind0 = h->values_to_mat[HEAP_MEM * i];
	i_ind1 = h->values_to_mat[HEAP_MEM * i + 1];

	tmp = h->mat_to_values[j_ind0][j_ind1];
	h->mat_to_values[j_ind0][j_ind1] = h->mat_to_values[i_ind0][i_ind1];
	h->mat_to_values[i_ind0][i_ind1] = tmp;

	set_key(h, j, i_k);
	set_indexes(h, j, i_ind0, i_ind1);
	set_key(h, i, j_k);
	set_indexes(h, i, j_ind0, j_ind1);
}


void heapify_down(matrice_max_heap *h, int i) {
	int l, r, max;
	l = left_loc(i);
	r = right_loc(i);
	max = i;
	if (l<heap_len(h) && left(h, i)>key(h, max))
		max = l;
	if (r<heap_len(h) && right(h, i)>key(h, max))
		max = r;
	if (max > i) {
		interswitch(h, i, max);
		heapify_down(h, max);
	}
}


void heapify_up(matrice_max_heap *h, int i) {
	while (i > 1 && key(h, i) > parent(h, i)) {
		interswitch(h, i, parent_loc(i));
		i = parent_loc(i);
	}
}


void heapify(matrice_max_heap *h) {
	int i, l;
	l = heap_len(h);
	for (i = l - 1; i >= 1; --i) {
		heapify_down(h, i);
	}
}

/*[i, j] where max value obtained in a_ij*/
void heap_max(matrice_max_heap *h, int arr[2]) {
	int i, j;
	i = h->values_to_mat[1 * HEAP_MEM];
	j = h->values_to_mat[1 * HEAP_MEM + 1];
	arr[0] = i;
	arr[1] = j;
}

void update_key(matrice_max_heap *h, int _ind, double v) {
	if (fabs(v) > key(h, _ind)) {
		set_key(h, _ind, fabs(v));
		heapify_up(h, _ind);
	}
	else if (fabs(v) < key(h, _ind)) {
		set_key(h, _ind, fabs(v));
		heapify_down(h, _ind);
	}
}

int init_max_heap_from_matrice(matrice_max_heap *h, double **mat,
	int **heap_locations_map, int n) {
	int i, j, k, m;
	m = (int)pow(n, 2) - n + 1;
	h->values = (double*)malloc(m * sizeof(double));
	if (h->values == NULL) return -1;
	h->values_to_mat = (int*)malloc(m * HEAP_MEM * sizeof(int));
	if (h->values_to_mat == NULL) return -1;

	h->mat_to_values = heap_locations_map;
	h->values[0] = -1;
	h->values_to_mat[0] = m;
	h->values_to_mat[1] = -1;
	k = 1;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			if (i == j)continue;
			h->values[k] = fabs(mat[i][j]);
			h->values_to_mat[k * HEAP_MEM] = i;
			h->values_to_mat[k * HEAP_MEM + 1] = j;
			heap_locations_map[i][j] = k;
			++k;
		}
	}
	heapify(h);


	return 0;
}

void free_heap(matrice_max_heap *h, int mat_dim) {
	int i;
	free(h->values);
	free(h->values_to_mat);
	for (i = 0; i < mat_dim; ++i) free(h->mat_to_values[i]);
	free(h->mat_to_values);
}
