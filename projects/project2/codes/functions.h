/*
Header file to declare functions.
*/

using namespace std;
using namespace arma;

void Find_MaxElement_and_MaxIndices(int&, int&, double&, mat, int);
void Compute_Trigonometric_Functions(int, int, mat, int, double&, double&, double&, double&);
mat FillUnitaryMatrix(int, int, int, double, double);
string OrthonormalityPreservationTest(mat, mat, int);
string ConservationOfEigenvalues(mat, vec, int);
