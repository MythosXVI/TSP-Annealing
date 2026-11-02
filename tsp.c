#include<limits.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<time.h>

#include "tsp.h"

#define MAX 1024
#define PI 3.141592
#define RRR 6378.388

int dimensions;

#define INX(d, i, j) ((d)[ (i) * (dimensions) + (j)])


/*---------------------------file parsing start----------------------------------*/
void read_tsp(char *file_name, tsp *t)
{
    int i;
    char line[MAX], *cp; /* cp is a pointer to work with */
    FILE *ifp = fopen(file_name, "r");
    t->dimensions = -1;
    t->type = NULL;

    getlin(ifp, line, MAX);

    /* Parse the header and get the relavent info */
    /* Use that info to initialize type tps */
    while(t->dimensions < 0 || t->type == NULL)
    {
        getlin(ifp, line, MAX);
        if(strstr(line, "DIMENSION"))
        {
            cp = strstr(line, ":");
            cp++;
            t->dimensions = atoi(cp);
        }else if(strstr(line, "EDGE_WEIGHT_TYPE"))
        {
            cp = strstr(line, ":");
            while(*++cp == ' ')
                ;
            t->type = malloc(strlen(cp) + 1);
            strcpy(t->type, cp);
        }
    }

    /* Allocate just enough space, now tsp knows its size as well */
    t->x = (double *)malloc( t->dimensions * sizeof(double));
    t->y = (double *)malloc( t->dimensions * sizeof(double));

    dimensions = t->dimensions;

    /* Skip other fields: name, etc. */
    while(strstr(line, "NODE_COORD_SECTION") == NULL)
        getlin(ifp, line, MAX);

    /* loop over each line, sscanf() is goated */
    i = 0;
    while(getlin(ifp, line, MAX) > 0)
    {
        sscanf(line, "%*d %lf %lf", &t->x[i], &t->y[i]);
        i++;
    }
    fclose(ifp);
}

int getlin(FILE *ifp, char s[], int max)
{
    int c, i;
    for(i = 0; i < max && (c = getc(ifp)) && c != '\n' && c != EOF; i++)
        s[i] = c;

    if(c == '\n'){
        s[i] = '\n';
        i++;
    }
    s[i] = '\0';
    return i;
}

/* ---------------------------Data Conversion---------------------------------- */


int nint(double num){return (int) (num + 0.5);};

int* EUC_2D(tsp* ts)
{
    int i,j;
    double xd, yd;
    int *matrix = (int *)malloc((ts->dimensions * ts->dimensions) * sizeof(int));

    for( i = 0; i < ts->dimensions; i++)
    {
        for(j = 0; j < ts->dimensions; j++){
            xd = ts->x[i] - ts->x[j];
            yd = ts->y[i] - ts->y[j];
            INX(matrix, i, j) = nint(sqrt(xd * xd + yd * yd));
        }
    }
    return matrix;
}

int* GEO(tsp *ts)
{
    int *matrix, i, j, deg;
    double latitude[ts->dimensions], longitude[ts->dimensions], min, q1, q2, q3;

    matrix = (int *) malloc((ts->dimensions * ts->dimensions) * sizeof(int));

    /* convert input to geographical latitude and longitude given in radians */
    for(i = 0; i < ts->dimensions; i++)
    {
        deg = (int)(ts->x[i]); /* might need to swich 'nint' for regular casting: (int) */
        min = ts->x[i] - deg;
        latitude[i] = PI * (deg + 5.0 * min / 3.0) /180;
        deg = (int)ts->y[i];
        min = ts->y[i] - deg;
        longitude[i] = PI * (deg + 5.0 * min / 3.0) /180;
    }

    for(i = 0; i < ts->dimensions; i++)
        for(j = 0; j < ts->dimensions; j++)
        {
            q1 = cos(longitude[i] - longitude[j]);
            q2 = cos(latitude[i] - latitude[j]);
            q3 = cos(latitude[i] + latitude[j]);
            INX(matrix, i, j) = (int) (RRR * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3) ) + 1.0);
        }

    return matrix;
}
/* ---------------------------Data Processing---------------------------------- */

int rand_point(int dimensions)
{
    return (int)(dimensions * (rand() / (RAND_MAX + 1.0)));
}

int not_visited(int *tour, int node, int filled)
{
    int i;
    for(i = 0; i < filled; i++)
    {
        if(node == tour[i])
            return 0;
    }
    return 1;
}


int *nearest_neighbor(int *matrix, int dimension, int start_city)
{
    int *tour, i, j, curr_city, min, closest = 0;
    if((tour = (int *)calloc(dimensions , sizeof(int))) == NULL)
    {
        printf("Error: nearest neighbor can't allocate error");
        return NULL;
    }

    tour[0] = start_city;

    for(i = 1; i < dimension; i++)
    {
        min = INT_MAX;
        for(j = 0; j < dimension; j++)
        {
            if(not_visited(tour, j, i))
            {
                curr_city = INX(matrix, tour[i - 1], j);
                if(curr_city < min){
                    min = curr_city;
                    closest = j;
                }
            }
        }
        tour[i] = closest;
    }

    return tour;
}

int tour_length(int *tour, int *matrix)
{
    int i, length = 0;
    for(i = 0; i < dimensions; i++)
        length += INX(matrix, tour[i], tour[(i + 1) % dimensions]);

    return length;
}

void two_opt(int *matrix, int *tour, int node1, int node2)
{
    int i, j;
    int mod_tour[dimensions];

    for(i = 0; i < dimensions; i++)
        mod_tour[i] = tour[i];

    for(j = node1,i = node2; i >= node1; i--, j++)
        mod_tour[j] = tour[i];

    for(i = 0; i < dimensions; i++)
        tour[i] = mod_tour[i];
}

int two_opt_delta(int *matrix, int *tour, int node1, int node2)
{
    int a, b, c, d, delta;
    a = tour[node1 - 1];
    b = tour[node1];
    c = tour[node2];
    d = tour[(node2 + 1) % dimensions];

    delta = INX(matrix, a, c) + INX(matrix, b, d) - INX(matrix, a, b) - INX(matrix, c, d);

    return delta;
};
