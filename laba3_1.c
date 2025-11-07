#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

double func(double x) {
  return 4 / (1 + x * x);
}

int main(int argc, char **argv) {
  int commsize, rank;
  const double a = 0.0;
  const double b = 1.0;
  const int n = 100;
  double h = (b - a) / n;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int points = n / commsize;
  int lb = rank * points;
  int ub = (rank == commsize - 1) ? (n - 1) : (lb + points - 1);

  double sum = 0.0;
  for (int i = lb; i <= ub; i++) {
    sum += func(a + h * (i + 0.5));
  }

  double gsum = 0.0;
  MPI_Reduce(&sum, &gsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    gsum *= h;
    printf("Result PI: %.12f\n", gsum);
  }

  MPI_Finalize();  
  return 0;
}