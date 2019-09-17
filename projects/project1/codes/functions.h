//Declaration of functions in header file.

void LU_thomas(double*, double*, double*, double*, double*, double*, double*, double*, double*, int);
void f(double, double, double&);
void LU_decomposition(double*, double*, double*, double*, double*, double*, int);
void Forward_substitutionLU(double*, double*, double*, int);
void Back_substitutionLU(double*, double*, double*, double*, int);
void Forward_substitution(double*, double*, double*, double*, int);
void Back_substitution(double*, double*, double*, double*, int);
void closed_form_solution(double, double&);
void compute_errors(double*, double*, double*, int);
void SpecialThomas(double*, double*, int);
