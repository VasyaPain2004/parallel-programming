#include <iostream>
#include <omp.h>

int main() {
    printf("=== ONLY PARALLEL (WRONG) ===\n");
    #pragma omp parallel
    {
        for (int i = 0; i < 5; ++i) {
            printf("Thread %d executes i=%d\n", omp_get_thread_num(), i);
        }
    }
    
    printf("\n=== PARALLEL FOR (CORRECT) ===\n");
    #pragma omp parallel for
    for (int i = 0; i < 5; ++i) {
        printf("Thread %d executes i=%d\n", omp_get_thread_num(), i);
    }
    
    return 0;
}