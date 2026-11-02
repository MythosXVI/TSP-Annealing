#include<string.h>
#include<stdio.h>

typedef struct _tsp{
	int dimensions;
	char *type;
	double *x;
	double *y;
}tsp;

extern int dimensions;

int getlin(FILE *, char [], int);
void read_tsp(char *file_name, tsp *);
int nint(double);
int* GEO(tsp *);
int* EUC_2D(tsp *);
int* nearest_neighbor(int *, int, int);
int rand_point(int);
int tour_length(int*, int*);
void two_opt(int *, int*, int node1, int node2);
int two_opt_delta(int*, int*, int, int);
