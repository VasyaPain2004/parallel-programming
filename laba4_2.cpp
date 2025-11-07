#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <cmath>

int main() {
    const long long TOTAL_POINTS = 1000000000;
    const double R = 5;
    const int NUM_THREADS = 8;
    
    double* thread_pi_values = new double[NUM_THREADS];
    long long* thread_points = new long long[NUM_THREADS];
    
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        int thread_id = omp_get_thread_num();
        long long local_points_in_circle = 0;
        
        unsigned int seed = time(NULL) + thread_id;
        
        long long points_per_thread = TOTAL_POINTS;
        
        for (long long i = 0; i < points_per_thread; i++) {
            double x = (double)rand_r(&seed) / RAND_MAX * R;
            double y = (double)rand_r(&seed) / RAND_MAX * R;
            
            if (x * x + y * y <= R * R) {
                local_points_in_circle++;
            }
        }
        
        double local_pi = 4.0 * (double)local_points_in_circle / points_per_thread;
        thread_pi_values[thread_id] = local_pi;
        thread_points[thread_id] = points_per_thread;
        
        printf("Thread %d: PI = %.10f, Points in circle = %lld/%lld\n", 
               thread_id, local_pi, local_points_in_circle, points_per_thread);
    }
    
    double sum_pi = 0.0;
    for (int i = 0; i < NUM_THREADS; i++) {
        sum_pi += thread_pi_values[i];
    }
    double average_pi = sum_pi / NUM_THREADS;
    
    printf("\n");
    for (int i = 0; i < NUM_THREADS; i++) {
        double error = fabs(thread_pi_values[i] - M_PI);
        printf("  Thread %d: PI = %.10f, Error = %.10f\n", 
               i, thread_pi_values[i], error);
    }
    printf("\n");
    printf("Average PI:    %.10f\n", average_pi);
    printf("Real PI:     %.10f\n", M_PI);
    printf("Average Error: %.10f\n", fabs(average_pi - M_PI));
    
    delete[] thread_pi_values;
    delete[] thread_points;
    
    return 0;
}