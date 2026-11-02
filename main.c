#include<stdio.h>
#include<string.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
#include <sys/time.h>

#include "tsp.h"

#define TRUE 1
#define FALSE 0

double accept_res(double delta, double temp);
double rand_num();


int main(int argc, char *argv[])
{
    int i, node1, node2, iter, accept, curr_length, best_length, tempo, optimal;
    int *matrix, *curr_tour;
    int *best_tour;
    double delt, temp, cooling, run_time;
    tsp ts;

    struct timeval start, end;
    gettimeofday(&start, NULL);

    srand((unsigned) time(NULL));
    if(argc == 3)
    {
        fprintf(stderr, "usage: tsp filename iteraions");
    }
    read_tsp(*++argv, &ts);

    if(strstr(ts.type, "EUC_2D"))
        matrix = EUC_2D(&ts);
    else if(strstr(ts.type, "GEO"))
        matrix = GEO(&ts);

    if(strstr(*argv, "burma"))
    {
        optimal = 3323;
    }else if(strstr(*argv, "berlin"))
    {
        optimal = 7542;
    }

    iter = atoi(*++argv);
    temp = 100;
    cooling = 0.999999;

    curr_tour = (int *) malloc(dimensions * sizeof(int) );
    best_tour = (int *) malloc(dimensions * sizeof(int) );
    curr_tour = nearest_neighbor(matrix, dimensions, 0);

    memcpy(best_tour, curr_tour, dimensions * sizeof(int));


    printf("\n\n-----------------Simulated Annealing----------------------\n");


    curr_length = best_length = tour_length(curr_tour, matrix);
    while(iter-- > 0)
    {
        do{
            node1 = rand_point(dimensions);
            node2 = rand_point(dimensions);

            if(node1 > node2){ tempo = node1; node1 = node2; node2 = tempo;}
        }while(node1 < 1 || node1 == node2);


        delt = two_opt_delta(matrix, curr_tour, node1, node2);
        if(delt < 0)
            accept = TRUE;
        else
            accept = rand_num() < accept_res(delt, temp);

        if(accept)
        {
            two_opt(matrix, curr_tour, node1, node2);
            curr_length += delt;
            if(curr_length < best_length)
            {
                best_length = curr_length;
                memcpy(best_tour, curr_tour, dimensions * sizeof(int));
            }
        }
        printf("\ni: %d | length: %d | temp: %.3f | delt: %.3f\n", iter, curr_length, temp, delt);
        temp *= cooling;
        if(temp < 1e-14)
            temp = 1e-14;
    }

    best_length = tour_length(best_tour, matrix);

    gettimeofday(&end, NULL);

    run_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;

    printf("-------------result---------------------------");
    printf("\nlength: %d | optimal: %d | diff: %d| runtime: %.4f \ntour: ", best_length, optimal, best_length - optimal, run_time);
    printf("[");
    for(i = 0; i < dimensions; i++)
        printf("%d ", best_tour[i]);
    printf("]\n");
}



double accept_res(double delta, double temp)
{
    return exp(-delta/temp);
}

double rand_num()
{
    return (double)rand() / RAND_MAX;
}
