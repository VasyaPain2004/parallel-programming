#include <iostream>
#include <omp.h>
using namespace std;

int main() {
    int N = 3;

    double** a = new double*[N];
    double* b = new double[N];
    double* x = new double[N];

    for (int i = 0; i < N; i++) {
        a[i] = new double[N];
    }

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            a[i][j] = (i + 1) * (j + 1) + 2.0 + (i == j ? N * 100.0 : 0);
        }
        b[i] = (i + 1) * 10.0; 
    }

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        cout << a[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl;

    for (int i = 0; i < N; i++) {
      cout << b[i] << " ";
    }
    cout << endl;

    double s;
    double t = omp_get_wtime();

    #pragma omp parallel
    {
      for (int k = 0; k < N - 1; k++) {
        double pivot = a[k][k];

        #pragma omp for schedule(dynamic)
        for (int i = k + 1; i < N; i++) {
          double lik = a[i][k] / pivot;
          for (int j = k; j < N; j++)
            a[i][j] -= lik * a[k][j];
          b[i] -= lik * b[k];
        }
      }
      for (int k = N - 1; k >= 0; k--) {
        s = 0;
        #pragma omp barrier

        #pragma omp for reduction(+:s)
        for (int i = k + 1; i < N; i++) {
          s += a[k][i] * x[i];
        }
        #pragma omp single 
        x[k] = (b[k] - s) / a[k][k];
      }
    }

    t = omp_get_wtime() - t;
    cout << "Time: " << t << endl;

    for (int i = 0; i < N; i++) {
      cout << x[i] << " ";
    }
    cout << endl;

    for (int i = 0; i < N; i++) {
        delete[] a[i];
    }

    delete[] a;
    delete[] b;
    delete[] x;
    
    return 0;
}