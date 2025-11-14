#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int N = 1000;
    const double epsilon = 1e-10; 
    const long int max_iter = 10000000000;

    if (N % size != 0) {
        if (rank == 0) printf("N must be divisible by number of processes\n");
        MPI_Finalize();
        return 1;
    }
    int local_n = N / size; 

    double* u_old = (double*)malloc((local_n + 2) * sizeof(double));
    double* u_new = (double*)malloc((local_n + 2) * sizeof(double));

    double h = 1.0 / N;
    for (int i = 0; i < local_n + 2; i++) {
        u_old[i] = 0.0;
    }

    if (rank == 0) {
        u_old[0] = 1.0;  
    }
    if (rank == size - 1) {
        u_old[local_n + 1] = 0.0; 
    }

    int iter = 0;
    double error = 1.0;

    double start_time = MPI_Wtime();

    while (error > epsilon && iter < max_iter) {
        MPI_Request requests[4];
        int req_count = 0;

        if (rank < size - 1) {
            MPI_Isend(&u_old[local_n], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &requests[req_count++]);
            MPI_Irecv(&u_old[local_n + 1], 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &requests[req_count++]);
        }

        if (rank > 0) {
            MPI_Isend(&u_old[1], 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &requests[req_count++]);
            MPI_Irecv(&u_old[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &requests[req_count++]);
        }

        MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);

        double local_error = 0.0;
        for (int i = 1; i <= local_n; i++) {
            u_new[i] = (u_old[i - 1] + u_old[i + 1]) / 2.0;
            double diff = fabs(u_new[i] - u_old[i]);
            if (diff > local_error) local_error = diff;
        }

        for (int i = 1; i <= local_n; i++) {
            u_old[i] = u_new[i];
        }

        MPI_Allreduce(&local_error, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        iter++;
    }

    double end_time = MPI_Wtime();

    double* u_global = NULL;
    if (rank == 0) {
        u_global = (double*)malloc((N + 1) * sizeof(double));
        u_global[0] = 1.0;
    }

    MPI_Gather(&u_old[1], local_n, MPI_DOUBLE, 
               rank == 0 ? &u_global[1] : NULL, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == size - 1) {
        if (rank == 0) u_global[N] = 0.0;
    } else if (rank == 0) {
        u_global[N] = 0.0; 
    }

    if (rank == 0) {
        for(int i=0; i<=N; i++) {
          printf("u[%d] = %f\n", i, u_global[i]);
        }
        printf("Time %f\n", end_time - start_time);
        printf("iter count %d\n", iter);
        free(u_global);
    }

    free(u_old);
    free(u_new);
    MPI_Finalize();
    return 0;
}