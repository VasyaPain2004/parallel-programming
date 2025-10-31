#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <cmath>

int main() {
    const long long TOTAL_POINTS = 10000000;
    const double R = 5;
    const int NUM_THREADS = 4;
    
    double start = omp_get_wtime();
    
    // Массивы для хранения результатов каждого потока
    double* thread_pi_values = new double[NUM_THREADS];
    long long* thread_points_in_circle = new long long[NUM_THREADS];
    long long* thread_total_points = new long long[NUM_THREADS];
    
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        int thread_id = omp_get_thread_num();
        long long local_points_in_circle = 0;
        
        unsigned int seed = time(NULL) + thread_id;
        long long points_per_thread = TOTAL_POINTS / NUM_THREADS;
        
        for (long long i = 0; i < points_per_thread; i++) {
            double x = (double)rand_r(&seed) / RAND_MAX * R;
            double y = (double)rand_r(&seed) / RAND_MAX * R;
            
            if (x * x + y * y <= R * R) {
                local_points_in_circle++;
            }
        }
        
        // Сохраняем сырые данные, а не вычисленное π
        thread_points_in_circle[thread_id] = local_points_in_circle;
        thread_total_points[thread_id] = points_per_thread;
        
        // Вычисляем π для этого потока (только для вывода)
        double local_pi = 4.0 * (double)local_points_in_circle / points_per_thread;
        thread_pi_values[thread_id] = local_pi;
        
        printf("Thread %d: PI = %.10f, Points in circle = %lld/%lld\n", 
               thread_id, local_pi, local_points_in_circle, points_per_thread);
    }
    
    // ПРАВИЛЬНОЕ вычисление среднего - через объединение всех точек
    long long total_points_in_circle = 0;
    long long total_actual_points = 0;
    
    for (int i = 0; i < NUM_THREADS; i++) {
        total_points_in_circle += thread_points_in_circle[i];
        total_actual_points += thread_total_points[i];
    }
    
    double correct_average_pi = 4.0 * (double)total_points_in_circle / total_actual_points;
    
    // Неправильное среднее арифметическое (для сравнения)
    double arithmetic_average_pi = 0.0;
    for (int i = 0; i < NUM_THREADS; i++) {
        arithmetic_average_pi += thread_pi_values[i];
    }
    arithmetic_average_pi /= NUM_THREADS;
    
    double end = omp_get_wtime();
    
    printf("\nIndividual PI values:\n");
    for (int i = 0; i < NUM_THREADS; i++) {
        double error = fabs(thread_pi_values[i] - M_PI);
        printf("  Thread %d: PI = %.10f, Error = %.10f\n", 
               i, thread_pi_values[i], error);
    }
    
    printf("\nComparison of averaging methods:\n");
    printf("================================\n");
    printf("Arithmetic average PI:  %.10f (Error = %.10f)\n", 
           arithmetic_average_pi, fabs(arithmetic_average_pi - M_PI));
    printf("Correct average PI:     %.10f (Error = %.10f)\n", 
           correct_average_pi, fabs(correct_average_pi - M_PI));
    printf("Actual PI:              %.10f\n", M_PI);
    printf("\nTotal points: %lld, Total in circle: %lld\n", 
           total_actual_points, total_points_in_circle);
    printf("Time: %.6f seconds\n", end - start);
    
    delete[] thread_pi_values;
    delete[] thread_points_in_circle;
    delete[] thread_total_points;
    
    return 0;
}