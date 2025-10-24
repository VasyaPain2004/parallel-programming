#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char** argv)
{
    int rank, size;
    int i, j, n = 3200;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n_partial = n / size;
    
    double *a_partial = malloc(n_partial * n * sizeof(double)); 
    double *x = malloc(n * sizeof(double)); 
    double *y_partial = malloc(n_partial * sizeof(double)); 
    double *y_total = malloc(n * sizeof(double)); 
    double *a = malloc(n * n * sizeof(double)); 
    
    if (a_partial == NULL || x == NULL || y_partial == NULL || y_total == NULL || a == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    if (rank == 0)
    {
        printf("Matrix\n");
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (i == j)
                    a[i*n + j] = 1;
                else
                    a[i*n + j] = 2;
                printf("%f\t", a[i*n + j]);
            }
            printf("\n");
        }
        printf("Vector\n");
        for (i = 0; i < n; i++)
        {
            x[i] = i + 1;
            printf("%f\t", x[i]);
        }
        printf("\n");
    }
    
    double t = MPI_Wtime();
    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
    MPI_Scatter(a, n_partial * n, MPI_DOUBLE, a_partial, n_partial * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for (i = 0; i < n_partial; i++)
    {
        y_partial[i] = 0.0;
        for (j = 0; j < n; j++)
            y_partial[i] += a_partial[i*n + j] * x[j];
    }
    
    MPI_Gather(y_partial, n_partial, MPI_DOUBLE, y_total, n_partial, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    t = MPI_Wtime() - t;
    
    if (rank == 0)
    {
        printf("time = %f\n", t);
        printf("Result vector:\n");
        for (i = 0; i < n; i++) {
            printf("%10.5f ", y_total[i]);
        }
        printf("\n");
    }
    
    free(a_partial);
    free(a);
    free(x);
    free(y_partial);
    free(y_total);

    MPI_Finalize();
    return 0;
}