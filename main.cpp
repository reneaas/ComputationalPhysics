#include <iostream>
#include <cmath>


using namespace std;

void f(double, double*);

int main(){
    double** A;
    int n = 10;
    A = new double*[n];
    for (int i=0; i < n; i++){
        A[i] = new double[n];
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            A[i][j] = 1.0;
        }
    }

    for (int i = 0; i < n; i++){
        delete[] A[i];
    }
    delete[] A;
    return 0;
}


