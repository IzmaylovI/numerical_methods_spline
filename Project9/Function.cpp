#include "Function.h"

double* Tridiagonal_matrix_algorithm(int n, double* A, double* B, double* C, double* D)  //A, B(n-1); C, D(n);  
{
    double* X = new double[n];
    double* alpha = new double[n - 1];
    double* beta = new double[n - 1];

    alpha[0] = B[0] / C[0];
    beta[0] = D[0] / C[0];

    for (int i = 1; i < n - 1; i++)
    {
        double temp = C[i] - A[i - 1] * alpha[i - 1];
        alpha[i] = B[i] / temp;
        beta[i] = (D[i] - A[i - 1] * beta[i - 1]) / temp;
    }

    X[n - 1] = (D[n - 1] - A[n - 2] * beta[n - 2]) / (C[n - 1] - A[n - 2] * alpha[n - 2]);

    for (int i = n - 2; i >= 0; i--)
    {
        X[i] = beta[i] - alpha[i] * X[i + 1];
    }

    delete[] alpha;
    delete[] beta;

    return X;
}

void interp_spline(int grid_dim, double step, double* X, double* Y, double* a, double* b, double* c, double* d) {

    double* A = new double[grid_dim - 1];
    double* B = new double[grid_dim - 1];
    double* C = new double[grid_dim];
    double* D = new double[grid_dim];


    // coeff Ai
    for (int i = 0; i < grid_dim - 2; i++) {
        A[i] = step;
    }
    A[grid_dim - 2] = 0;

    // coeff Ci
    C[0] = 1;
    for (int i = 1; i < grid_dim - 1; i++) {
        C[i] = 2 * (step + step);
    }
    C[grid_dim - 1] = 1;

    //coeff Bi
    B[0] = 0;
    for (int i = 1; i < grid_dim - 1; i++) {
        B[i] = step;
    }

    //coeff Di
    D[0] = 0;
    for (int i = 1; i < grid_dim - 1; i++) {
        D[i] = 6 * ((Y[i + 1] - Y[i]) / step - (Y[i] - Y[i - 1]) / step);
    }
    D[grid_dim - 1] = 0;

    //coeff ci
    double* c_tmp = Tridiagonal_matrix_algorithm(grid_dim, A, B, C, D);

    for (int i = 0; i < grid_dim; i++) {
        c[i] = c_tmp[i];
    }
    delete[] c_tmp;

    // coeff ai
    for (int i = 0; i < grid_dim; ++i) {
        a[i] = Y[i + 1];
    }

    //coeff di
    for (int i = 0; i < grid_dim - 1; i++) {
        d[i] = (c[i + 1] - c[i]) / step;
    }

    //coeff bi
    for (int i = 0; i < grid_dim - 1; i++) {
        b[i] = (Y[i + 1] - Y[i]) / step + c[i + 1] * step / 3 + c[i] * step / 6;
    }

    delete[] A;
    delete[] B;
    delete[] C;
}