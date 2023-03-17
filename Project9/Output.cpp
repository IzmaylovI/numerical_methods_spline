#include "Output.h"
#include <iostream>

void Table1(std::ofstream& out, int grid_dim, double* X, 
    double* a, double* b, double* c, double* d) {

    //Table1 << "i, X(i-1), X(i), a(i), b(i), c(i), d(i)\n";

    for (int i = 0; i < grid_dim - 1; i++) {
        out << i + 1 << ',' << X[i] << ',' << X[i + 1] << ',' << a[i] << ',' << b[i] << ',' << c[i + 1] << ',' << d[i] << "\n";
    }
}

void Table2(std::ofstream& out, int N, double* X_control, double* Y_control, 
    double* Y_control_pr1, double* Y_control_pr2, double* S_control, 
    double* S_control_pr1, double* S_control_pr2) {
    
    //Table2 << "j,Xj,F(xj),S(xj),F(xj)-S(xj),F'(xj),S'(xj),F'(xj)-S'(xj),F''(xj),S''(xj),F''(xj)-S''(xj)\n";

    for (int i = 0; i < N; i++) {
        out << i << ',' << X_control[i] << ',' << Y_control[i] << ',' << S_control[i] << ',' << Y_control[i] - S_control[i] << ','
            << Y_control_pr1[i] << ',' << S_control_pr1[i] << ',' << Y_control_pr1[i] - S_control_pr1[i] << '\n';
    }
}

void Directory(std::ofstream& out, double grid_dim, int N, double S_Y_error, double S_Y_pr1_error, double S_Y_pr2_error) {
    
    out << grid_dim << ',' << N << ',' << S_Y_error << ','
        << S_Y_pr1_error << ',' << S_Y_pr2_error;
}