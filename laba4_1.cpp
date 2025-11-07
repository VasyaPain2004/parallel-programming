#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>

int main(int argc, char** argv) {
    long long n = 10000000;
    double sum = 0.0;
    
    #pragma omp parallel for reduction(+:sum)
    for (long long i = 0; i < n; ++i) {
        double denominator = 2.0 * i + 1.0;
        double term = 1.0 / denominator;
        
        if (i % 2 == 1) {
            term = -term;
        }
        
        sum += term;
    }
    
    double pi = 4 * sum;
    
    std::cout << std::setprecision(15);
    std::cout << "PI: " << pi << std::endl;
    std::cout << "Real PI: " << 3.14159265358979 << std::endl;
    // g++ -Wall -fopenmp example.cpp -o example
    // OMP_NUM_THREADS=3 ./laba4_1
    return 0;
}