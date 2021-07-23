double*** calcJacobi(double **, int);
double ** buildPij(double **, int, int, int, double, double);
int* getPivotIndexes(double **, int);
double* calcCandS(double **, int, int);
double ** MatriceMultiplication(double **, double **, int);
double ** calcAtag(double **, int, int, int, double, double);
int isDiagonalEnough(double **, double **, int);