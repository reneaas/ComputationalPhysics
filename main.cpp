#include <iostream>
#include <cmath>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char * argv[]){
    mat A =  mat(3,3);
    A.fill(1);
    A(0,0) = 1;
    A(0,1) = 2;
    A(0,2) = -1;
    A(1,0) = 2;
    A(1,1) = 3;
    A(1,2) = -3;
    A(2,0) = -1;
    A(2,1) = 2;
    A(2,2) = 3;
    cout << A << endl;

    for (int m = 0; m < 3; m++){
        for (int j = m; j < 3; j++){
            for (int k = m; k < 3; k++){
                A(j,k) = A(j,k) - (A(j,m)*A(k,m))/A(m,m);
            }
        }
    }
    cout << A << endl;
    return 0;
}
