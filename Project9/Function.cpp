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

void interp_spline(int grid_dim, double step, double* X, double* Y, double* a, double* b, double* c, double* d, double mu1, double mu2) {

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
    D[0] = mu1;
    for (int i = 1; i < grid_dim - 1; i++) {
        D[i] = 6 * ((Y[i + 1] - Y[i]) / step - (Y[i] - Y[i - 1]) / step);
    }
    D[grid_dim - 1] = mu2;

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

void get_function(int var, int grid_dim, double& step, double* X, double* Y) {

    switch (var)
    {
    case 0: {

        double l = -1;
        double r = 1;

        step = (r - l) / (grid_dim - 1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = -1 + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            if (X[i] <= 0) {
                Y[i] = X[i] * X[i] * X[i] + 3 * X[i] * X[i];
            }
            else {
                Y[i] = -1 * X[i] * X[i] * X[i] + 3 * X[i] * X[i];
            }
        }
        break;
    }
    case 1: {

        double l = 2;
        double r = 4;

        step = (r - l) / (grid_dim - 1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = log(X[i] + 1) / X[i];
        }

        break;
    }
    case 2: {

        double l = 1;
        double r = 3.14159265358979323846;

        step = (r - l) / (grid_dim - 1);


        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }
        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = sin(X[i] + 1) / X[i];
        }

        break;
    }
    case 3: {

        double l = 1;
        double r = 3.14159265358979323846;

        step = (r - l) / (grid_dim - 1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = (sin(X[i]) * sin(X[i])) / X[i];
        }

        break;
    }
    case 4: {

        double l = 0;
        double r = 3.14159265358979323846;

        step = (r - l) / (grid_dim - 1);;

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = X[i] * sin(X[i]) / 3;
        }

        break;
    }
    case 5: {

        double l = 2;
        double r = 4;

        step = (r - l) / (grid_dim - 1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = log(X[i] + 1) / X[i] + cos(10 * X[i]);
        }

        break;
    }
    case 6: {

        double l = 1;
        double r = 3.14159265358979323846;

        step = (r - l) / (grid_dim - 1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = sin(X[i] + 1) / X[i] + cos(10 * X[i]);
        }

        break;
    }
    case 7: {

        double l = 1;
        double r = 3.14159265358979323846;

        step = (r - l) / (grid_dim - 1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = (sin(X[i]) * sin(X[i])) / X[i] + cos(10 * X[i]);
        }

        break;
    }
    case 8: {

        double l = 0;
        double r = 3.14159265358979323846;

        step = (r - l) / (grid_dim - 1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = X[i] * sin(X[i]) / 3 + cos(10 * X[i]);
        }

        break;
    }
    }
}

void spline_znachenia(int grid_dim, double* S, double* S2, double* X, double* a, double* b, double* c, double* d) {
    S[0] = a[0] + b[0] * (X[0] - X[1]) + c[1] / 2 * (X[0] - X[1]) * (X[0] - X[1]) + d[0] / 6 * (X[0] - X[1]) * (X[0] - X[1]) * (X[0] - X[1]);
    for (int i = 1; i < grid_dim; i++) {
        S[i] = a[i - 1];
    }

    S2[0] = b[1] * X[0] + c[1] * (X[0] - X[1]) + d[1] / 2 * (X[0] - X[1]) * (X[0] - X[1]);
    for (int i = 1; i < grid_dim; i++) {
        S2[i] = b[i] * X[i];
    }
}

void test_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2){
    for (int i = 0; i < N; i++) {
        X_control[i] = -1 + step2 * i;
    }

    for (int i = 0; i < N; i++) {
        if (X_control[i] <= 0)
            Y_control[i] = X_control[i] * X_control[i] * X_control[i] + 3 * X_control[i] * X_control[i];
        else
            Y_control[i] = -X_control[i] * X_control[i] * X_control[i] + 3 * X_control[i] * X_control[i];
    }

    // первая производная 
    for (int i = 0; i < N; i++) {
        if (X_control[i] <= 0)
            Y_control_pr1[i] = 3 * X_control[i] * X_control[i] + 6 * X_control[i];
        else
            Y_control_pr1[i] = -3 * X_control[i] * X_control[i] + 6 * X_control[i];
    }
    

    // вторая производная 
    for (int i = 0; i < N; i++) {
        if (X_control[i] <= 0)
            Y_control_pr2[i] = 6 * X_control[i] + 6;
        else
            Y_control_pr2[i] = -6 * X_control[i] + 6;
    }
}

void var1_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2) {
    for (int i = 0; i < N; i++)
        X_control[i] = 2 + step2 * i;

    for (int i = 0; i < N; i++) {
        Y_control[i] = log(X_control[i] + 1) / X_control[i];

        Y_control_pr1[i] = 1 / (X_control[i] + X_control[i] * X_control[i]) - log(X_control[i] + 1) / (X_control[i] * X_control[i]);

        Y_control_pr2[i] = (-1 - 2 * X_control[i]) / ((X_control[i] + X_control[i] * X_control[i]) * (X_control[i] + X_control[i] * X_control[i])) - 1 / (X_control[i] * X_control[i] * (X_control[i] + 1)) + 2 * log(X_control[i] + 1) / (X_control[i] * X_control[i] * X_control[i]);

    }
}

void var2_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2) {
    for (int i = 0; i < N; i++)
        X_control[i] = 1 + step2 * i;

    for (int i = 0; i < N; i++) {
        Y_control[i] = sin(X_control[i] + 1) / X_control[i];

        Y_control_pr1[i] = (X_control[i] * cos(X_control[i] + 1) - sin(X_control[i] + 1)) / (X_control[i] * X_control[i]);

        Y_control_pr2[i] = -((X_control[i] * X_control[i] - 2) * sin(X_control[i] + 1) + 2 * X_control[i] * cos(X_control[i] + 1)) / (X_control[i] * X_control[i] * X_control[i]);

    }
}

void var3_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2) {
    for (int i = 0; i < N; i++)
        X_control[i] = 1 + step2 * i;

    for (int i = 0; i < N; i++) {
        Y_control[i] = sin(X_control[i]) * sin(X_control[i]) / X_control[i];

        Y_control_pr1[i] = (sin(X_control[i] * 2) * X_control[i] - sin(X_control[i]) * sin(X_control[i]))/(X_control[i] * X_control[i]);

        Y_control_pr2[i] = -(2 * X_control[i] * sin(2 * X_control[i]) - 2 * X_control[i] * X_control[i] * cos(2 * X_control[i]) - 2 * sin(X_control[i]) * sin(X_control[i])) / (X_control[i] * X_control[i] * X_control[i]);
    }
}

void var4_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2) {
    for (int i = 0; i < N; i++)
        X_control[i] = 0 + step2 * i;

    for (int i = 0; i < N; i++) {
        Y_control[i] = X_control[i] * sin(X_control[i]) / 3;

        Y_control_pr1[i] = (sin(X_control[i]) + X_control[i] * cos(X_control[i])) / 3;

        Y_control_pr2[i] = (2 * cos(X_control[i]) - X_control[i] * sin(X_control[i])) / 3;
    }
}

void var5_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2) {
    for (int i = 0; i < N; i++)
        X_control[i] = 2 + step2 * i;

    for (int i = 0; i < N; i++) {
        Y_control[i] = log(X_control[i] + 1) / X_control[i] + cos(10 * X_control[i]);

        Y_control_pr1[i] = 1 / (X_control[i] + X_control[i] * X_control[i]) - log(X_control[i] + 1) / (X_control[i] * X_control[i]) - 10 * sin(10 * X_control[i]);

        Y_control_pr2[i] = (-1 - 2 * X_control[i]) / ((X_control[i] + X_control[i] * X_control[i]) * (X_control[i] + X_control[i] * X_control[i])) - 1 / (X_control[i] * X_control[i] * (X_control[i] + 1)) + 2 * log(X_control[i] + 1) / (X_control[i] * X_control[i] * X_control[i]) - 100 * cos(10 * X_control[i]);

    }
}

void var6_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2) {
    for (int i = 0; i < N; i++)
        X_control[i] = 1 + step2 * i;

    for (int i = 0; i < N; i++) {
        Y_control[i] = sin(X_control[i] + 1) / X_control[i] + cos(10 * X_control[i]);

        Y_control_pr1[i] = (X_control[i] * cos(X_control[i] + 1) - sin(X_control[i] + 1)) / (X_control[i] * X_control[i]) - 10 * sin(10 * X_control[i]);

        Y_control_pr2[i] = -((X_control[i] * X_control[i] - 2) * sin(X_control[i] + 1) + 2 * X_control[i] * cos(X_control[i] + 1)) / (X_control[i] * X_control[i] * X_control[i]) - 100 * cos(10 * X_control[i]);
    }
}

void var7_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2) {
    for (int i = 0; i < N; i++)
        X_control[i] = 1 + step2 * i;

    for (int i = 0; i < N; i++) {
        Y_control[i] = sin(X_control[i]) * sin(X_control[i]) / X_control[i] + cos(10 * X_control[i]);

        Y_control_pr1[i] = (sin(X_control[i] * 2) * X_control[i] - sin(X_control[i]) * sin(X_control[i])) / (X_control[i] * X_control[i]) - 10 * sin(10 * X_control[i]);

        Y_control_pr2[i] = -(2 * X_control[i] * sin(2 * X_control[i]) - 2 * X_control[i] * X_control[i] * cos(2 * X_control[i]) - 2 * sin(X_control[i]) * sin(X_control[i])) / (X_control[i] * X_control[i] * X_control[i]) - 100 * cos(10 * X_control[i]);
    }
}

void var8_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2) {
    for (int i = 0; i < N; i++)
        X_control[i] = 0 + step2 * i;

    for (int i = 0; i < N; i++) {
        Y_control[i] = X_control[i] * sin(X_control[i]) / 3 + cos(10 * X_control[i]);

        Y_control_pr1[i] = (sin(X_control[i]) + X_control[i] * cos(X_control[i])) / 3 - 10 * sin(10 * X_control[i]);

        Y_control_pr2[i] = (2 * cos(X_control[i]) - X_control[i] * sin(X_control[i])) / 3 - 100 * cos(10 * X_control[i]);
    }
}

void spline_pr(int grid_dim, int N, int mu1, int mu2, double* X_control, double* S_control, double* S_control_pr1, double* S_control_pr2, double* a, double* b, double* c, double* d) {
    
    S_control[0] = a[0] + b[0] * (X_control[0] - X_control[3]) + c[1] / 2 * (X_control[0] - X_control[3]) * (X_control[0] - X_control[3]) + d[0] / 6 * (X_control[0] - X_control[3]) * (X_control[0] - X_control[3]) * (X_control[0] - X_control[3]);
    for (int i = 0; i < grid_dim-1; i++) {
        for (int j = 1; j <= 3; j++) {
            S_control[i * 3 + j] = a[i] + b[i] * (X_control[i * 3 + j] - X_control[(i + 1) * 3]) + c[i + 1] / 2 * (X_control[i * 3 + j] - X_control[(i + 1) * 3]) * (X_control[i * 3 + j] - X_control[(i + 1) * 3]) + d[i] / 6 * (X_control[i * 3 + j] - X_control[(i + 1) * 3]) * (X_control[i * 3 + j] - X_control[(i + 1) * 3]) * (X_control[i * 3 + j] - X_control[(i + 1) * 3]);
        }
    }

    S_control_pr1[0] = b[0] + c[1] * (X_control[0] - X_control[3]) + d[0] / 2 * (X_control[0] - X_control[3]) * (X_control[0] - X_control[3]);
    for (int i = 0; i < grid_dim-1; i++) {
        for (int j = 1; j <= 3; j++)
            S_control_pr1[i * 3 + j] = b[i] + c[i + 1] * (X_control[i * 3 + j] - X_control[(i + 1) * 3]) + d[i] / 2 * (X_control[i * 3 + j] - X_control[(i + 1) * 3]) * (X_control[i * 3 + j] - X_control[(i + 1) * 3]);
    }

    S_control_pr2[0] = mu1;
    for (int i = 0; i < grid_dim-1; i++) {
        for (int j = 1; j <= 3; j++)
            S_control_pr2[i * 3 + j] = c[i + 1] + d[i] * (X_control[i * 3 + j] - X_control[(i + 1) * 3]);
    }
    S_control_pr2[N - 1] = mu2;
}

void error(int N, double* Y_control, double* Y_control_pr1, double* Y_control_pr2, double* S_control, double* S_control_pr1, double* S_control_pr2, double& S_Y_error, double& S_Y_pr1_error, double& S_Y_pr2_error, int& X1, int& X2, int& X3) {
    for (int i = 0; i < N; ++i) {
        if (abs(Y_control[i] - S_control[i]) > S_Y_error) {
            S_Y_error = abs(Y_control[i] - S_control[i]);
            X1 = i;
        }
        if (abs(Y_control_pr1[i] - S_control_pr1[i]) > S_Y_pr1_error) {
            S_Y_pr1_error = abs(Y_control_pr1[i] - S_control_pr1[i]);
            X2 = i;
        }
        if (abs(Y_control_pr2[i] - S_control_pr2[i] > S_Y_pr2_error)) {
            S_Y_pr2_error = abs(Y_control_pr2[i] - S_control_pr2[i]);
            X3 = i;
        }
    }
}