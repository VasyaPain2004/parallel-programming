#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int ProcNum, ProcRank;

void Flip(double *B, int size) {
    double *temp = (double*)malloc(size * size * sizeof(double));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            temp[j * size + i] = B[i * size + j];
        }
    }
    for (int i = 0; i < size * size; i++) {
        B[i] = temp[i];
    }
    free(temp);
}

void MatrixMultiplicationMPI(double *A, double *B, double *C, int size) {
    int dim = size;
    int i, j, k, p, ind;
    double temp;
    MPI_Status Status;
    int ProcPartSize = dim / ProcNum; 
    int ProcPartElem = ProcPartSize * dim; 
    double* bufA = (double*)malloc(ProcPartElem * sizeof(double));
    double* bufB = (double*)malloc(ProcPartElem * sizeof(double));
    double* bufC = (double*)malloc(ProcPartElem * sizeof(double));
    int ProcPart = dim / ProcNum, part = ProcPart * dim;
    
    if (ProcRank == 0) {
        Flip(B, size);
    }
    
    MPI_Scatter(A, part, MPI_DOUBLE, bufA, part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(B, part, MPI_DOUBLE, bufB, part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    temp = 0.0;
    for (i = 0; i < ProcPartSize; i++) {
        for (j = 0; j < ProcPartSize; j++) {
            for (k = 0; k < dim; k++) 
                temp += bufA[i*dim+k] * bufB[j*dim+k];
            bufC[i*dim+j+ProcPartSize*ProcRank] = temp;
            temp = 0.0;
        }
    }

    int NextProc; int PrevProc;
    for (p = 1; p < ProcNum; p++) {
        NextProc = ProcRank + 1;
        if (ProcRank == ProcNum - 1) 
            NextProc = 0;
        PrevProc = ProcRank - 1;
        if (ProcRank == 0) 
            PrevProc = ProcNum - 1;
        MPI_Sendrecv_replace(bufB, part, MPI_DOUBLE, NextProc, 0, PrevProc, 0, MPI_COMM_WORLD, &Status);
        temp = 0.0;
        for (i = 0; i < ProcPartSize; i++) {
            for (j = 0; j < ProcPartSize; j++) {
                for (k = 0; k < dim; k++) {
                    temp += bufA[i*dim+k] * bufB[j*dim+k];
                }
                if (ProcRank - p >= 0 ) 
                    ind = ProcRank - p;
                else ind = (ProcNum - p + ProcRank);
                bufC[i*dim+j+ind*ProcPartSize] = temp;
                temp = 0.0;
            }
        }
    }
    
    MPI_Gather(bufC, ProcPartElem, MPI_DOUBLE, C, ProcPartElem, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(bufA);
    free(bufB);
    free(bufC);
}


int main(int argc, char **argv) { 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    
    int n = 6;
    
    if (n % ProcNum != 0) {
        if (ProcRank == 0) {
            printf("Error: Matrix size (%d) must be divisible by number of processes (%d)\n", n, ProcNum);
        }
        MPI_Finalize();
        return 1;
    }
    
    double *A = NULL;
    double *B = NULL;
    double *C = NULL;
    
    if (ProcRank == 0) {
        A = (double*)malloc(n * n * sizeof(double));
        B = (double*)malloc(n * n * sizeof(double));
        C = (double*)malloc(n * n * sizeof(double));
        
        srand(time(NULL));
        
        printf("Matrix A:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i * n + j] = (double)rand() / RAND_MAX * 10.0;
                printf("%6.2f ", A[i * n + j]);
            }
            printf("\n");
        }
        
        printf("\nMatrix B:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                B[i * n + j] = (double)rand() / RAND_MAX * 10.0;
                printf("%6.2f ", B[i * n + j]);
            }
            printf("\n");
        }
    } else {
        C = (double*)malloc(n * n * sizeof(double));
    }
    
    double start_time = MPI_Wtime();
    MatrixMultiplicationMPI(A, B, C, n);
    double end_time = MPI_Wtime();

    if (ProcRank == 0) {
        printf("\nResult C:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%8.2f ", C[i * n + j]);
            }
            printf("\n");
        }
        printf("\n");
        printf("Time: %f", end_time - start_time);
        printf("\n");
        
        free(A);
        free(B);
    }
    free(C);
    
    MPI_Finalize();
    return 0;
}


// import numpy as np

// # Матрица A
// A = np.array([
//   [6.77,  5.63 , 1.84, 8.15, 6.70, 0.60],
//   [0.97,  3.12,  0.67, 0.85, 3.99, 2.80],
//   [1.57,  2.41 , 3.60, 9.94, 7.43, 9.55],
//   [0.26,  0.39,  6.70, 0.55, 7.07, 6.25],
//   [1.22,  7.62,  4.76, 5.28, 0.66, 6.21],
//   [2.22, 7.43,  1.84, 4.06, 5.58, 8.54],
// ])

// # Матрица B
// B = np.array([
//     [4.66  , 6.55 , 1.66,  5.33 , 7.40,  5.65],
//   [8.13  , 8.97 , 8.06 , 1.72 , 8.92 , 5.49],
//   [1.28 ,  9.18,  5.88 , 7.98,  9.73,  2.95],
//   [4.23  , 0.94 , 0.57 , 8.99 , 6.22,  1.22],
//   [5.20  , 8.44 , 8.65 , 7.04,  2.50 , 4.23],
//   [5.58  , 7.16 , 0.78 , 7.24 , 2.49 , 8.18],
// ])
// # Умножение матриц

// C = np.dot(A, B)

// print("Matrix A:")
// print(A)
// print("\nMatrix B:")
// print(B)
// print("\nResult A × B:")
// print(C)