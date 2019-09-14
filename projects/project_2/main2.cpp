#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstdlib>


using namespace std;



//Declaration of functions

double** CreateMatrix(int);
void DestroyMatrix(double**, int);
void MatrixMultiplication(double**, double**, double**, int);
void SquareAllElements(double**, double**, int);
void Find_MaxElement_and_MaxIndices(int&, int&, double&, double**, int);
void Compute_Trigonometric_Functions(int, int, double**, int, double&, double&, double&, double&);
void FillUnitaryMatrix(double**, double**, int, int, int, double, double);



//Main program

int main(int argc, char* argv[]){

    //Declaration of variables:
    int n, RowIndex, ColumnIndex, k, l;         //Integers
    double sinus, cosinus, tangens, tau, tolerance, h, a, d, max_element;      //Floating points.
    double **A, **AA, **S, **S_transpose, **Storage_matrix, **B;      //Matrices.

    //Specify integers:
    //n = atoi(argv[1]);
    //cin >> n; //Temporary solution to specify the number n from "terminal".
    n = 3;
    //Specify floats:
    h = 1.0/((double) n);
    d = 2.0/(h*h);
    a = -1.0/(h*h);
    tolerance = 1e-9;
    max_element = 1.0;    //initial value to pass first check in while loop.
    int max_iterations = 2;


    //Create matrices:
    A = CreateMatrix(n);
    B = CreateMatrix(n);
    AA = CreateMatrix(n);
    S = CreateMatrix(n);
    S_transpose = CreateMatrix(n);
    Storage_matrix = CreateMatrix(n);

    //----------------Test zone ahead!--------------------

    //Fill matrices with random values to test function:
    //Fill up the tridiagonal matrix A:
    for (int i = 0; i < n; i++){
            if (i < n-1){
                A[i][i] = d;
                A[i][i+1] = a;
                A[i+1][i] = a;
            }
            else{
                A[i][i] = d;
            }
        }
    cout << "Initial matrix A" << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << A[i][j] << " ";
        }
        cout << " " << endl;
    }
    int iterations = 0;
    while (max_element*max_element > tolerance && iterations < max_iterations){
        //SquareAllElements(A, AA, n);
        Find_MaxElement_and_MaxIndices(RowIndex, ColumnIndex, max_element, A, n);

        cout << "max_element = " << max_element << endl;
        cout << "Row_index = " << RowIndex << endl;
        cout << "ColumnIndex = " << ColumnIndex << endl;
        //Next step: Compute trigonometric functions.
        Compute_Trigonometric_Functions(RowIndex, ColumnIndex, A, n, tau, tangens, cosinus, sinus);
        //Next step after that: Compute the unitary matrix
        k = RowIndex;
        l = ColumnIndex;
        FillUnitaryMatrix(S, S_transpose, k, l, n, cosinus, sinus);

        //Then compute the similary matrix and implement the whole charade in a while loop.
        MatrixMultiplication(Storage_matrix, A, S, n);
        MatrixMultiplication(B, S_transpose, Storage_matrix, n);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                A[i][j] = B[i][j];
            }
        }
        DestroyMatrix(Storage_matrix, n);
        DestroyMatrix(B,n);
        Storage_matrix = CreateMatrix(n);
        B = CreateMatrix(n);
        cout << "Matrix A = " << endl;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                cout << A[i][j] << " ";
            }
            cout << " " << endl;
        }
        iterations += 1;
    }

    return 0;
}

//Function specifications belows:


double **CreateMatrix(int n){
    double **mat;
    mat = new double*[n];
    for (int i = 0; i < n; i++){
        mat[i] = new double[n];
    }
    return mat;
}

void DestroyMatrix(double **mat, int n){
    for (int i = 0; i < n; i++){
        delete[] mat[i];
    }
    delete[] mat;
}

void MatrixMultiplication(double **C, double **A, double **B, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            for (int k = 0; k < n; k++){
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}

void SquareAllElements(double **A, double **AA, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i != j){
                AA[i][j] = A[i][j]*A[i][j];
            }
            else{
                AA[i][j] = 0.0;
            }
        }
    }
}


void Find_MaxElement_and_MaxIndices(int &RowIndex, int &ColumnIndex, double &max_element, double **A, int n){
    max_element = 0.0;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (abs(A[i][j]) > max_element && i != j){
                RowIndex = i;
                ColumnIndex = j;
                max_element = abs(A[i][j]);
            }
        }
    }
}


void Compute_Trigonometric_Functions(int row_index, int column_index, double **A, int n, double &tau, double &tangens, double &cosinus, double &sinus){
    if (A[row_index][column_index] != 0.0) {
        tau = (A[column_index][column_index] - A[row_index][row_index]) / (2.0 * A[row_index][column_index]);
        if (tau >= 0) {
            tangens = 1.0/(tau + sqrt(1.0 + tau * tau));
        } else {
            tangens = -1.0/(-tau + sqrt(1.0 + tau * tau));
        }
        cosinus = 1.0 / sqrt(1.0 + tangens * tangens);
        sinus = tangens*cosinus;
    }
    else{
        cosinus = 1.0;
        sinus = 0.0;
    }
}

void FillUnitaryMatrix(double **S, double **S_transpose, int k, int l, int n, double cosinus, double sinus){
    //First make the matrix the identity matrix.
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j){
                S[i][j] = 1.0;
                S_transpose[i][j] = 1.0;
            }
            else{
                S[i][j] = 0.0;
                S_transpose[i][j] = 0.0;
            }
        }
    }
    //Fill in the rest of the elements to make it into a rotation matrix.
    S[k][k] = cosinus;
    S[l][l] = cosinus;
    S[k][l] = -sinus;
    S[l][k] = sinus;

    //Fill the transposed unitary matrix.
    S_transpose[k][k] = cosinus;
    S_transpose[l][l] = cosinus;
    S_transpose[l][k] = -sinus;
    S_transpose[k][l] = sinus;
}
