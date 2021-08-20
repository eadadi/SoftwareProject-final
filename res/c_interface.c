#include "value_vector_map_struct.h"
#define NORMAL_LENGTH 10
#define MAX_NUM_OF_POINTS 1000

#ifndef C_ENUM_GOAL
#define C_ENUM_GOAL

enum goal { wam, ddg, lnorm, jacobi, spk };
#endif
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
	int i, cnt;
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
int getdata(FILE *f, double ***data, int *features_num) {
	double *vector;
	char *line;
	int i, len;
	*data = (double**)malloc(MAX_NUM_OF_POINTS * sizeof(double*));
	len = freadline(&line, f);
	if (len == -1) return -1;
	*features_num = countCommas(line) + 1;
	i = 0;
	while (line[0] != '\0') {
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
