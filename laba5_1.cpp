#include <iostream>
#include <string>
#include <omp.h>
using namespace std;

bool checkPassword(const string& password) {
    return password == "zzzr";
}

int main() {
    omp_set_num_threads(4);

    string password;

    double start_time = omp_get_wtime();
    
    #pragma omp parallel for
    for (char c1 = 'a'; c1 <= 'z'; c1++) {
        for (char c2 = 'a'; c2 <= 'z'; c2++) {
            for (char c3 = 'a'; c3 <= 'z'; c3++) {
                for (char c4 = 'a'; c4 <= 'z'; c4++) {
                    string candidate = {c1, c2, c3, c4};
                    
                    if (checkPassword(candidate)) {
                        password = candidate;
                    }
                }
            }
        }
    }

    double end_time = omp_get_wtime();
    
    if (!password.empty()) {
        cout << "Password: " << password << endl;
    }

    cout << "Time: " << end_time - start_time << endl;
    
    return 0;
}