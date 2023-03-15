#include "Output.h"

void Table1(int grid_dim, double* X, double* a, double* b, double* c, double* d) {
    std::ofstream Table1("Table_1.csv");
    Table1 << "i, X(i-1), X(i), a(i), b(i), c(i), d(i)\n";

    for (int i = 0; i < grid_dim - 1; i++) {
        Table1 << i + 1 << ',' << X[i + 1] << ',' << X[i] << ',' << a[i] << ',' << b[i] << ',' << c[i + 1] << ',' << d[i] << "\n";
    }
}

void Table2(int N, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2, double* S_control, double* S_control_pr1, double* S_control_pr2) {
    std::ofstream Table2("Table_2.csv");
    
    Table2 << "j,Xj,F(xj),S(xj),F(xj)-S(xj),F'(xj),S'(xj),F'(xj)-S'(xj),F''(xj),S''(xj),F''(xj)-S''(xj)\n";

    for (int i = 0; i < N; i++) {
        Table2 << i << ',' << X_control[i] << ',' << Y_control[i] << ',' << S_control[i] << ',' << Y_control[i] - S_control[i] << ','
            << Y_control_pr1[i] << ',' << S_control_pr1[i] << ',' << Y_control_pr1[i] - S_control_pr1[i] << ','
            << Y_control_pr2[i] << ',' << S_control_pr2[i] << ',' << Y_control_pr2[i] - S_control_pr2[i] << '\n';
    }
}