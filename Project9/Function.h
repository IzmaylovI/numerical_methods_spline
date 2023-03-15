#pragma once
#include <cmath>


double* Tridiagonal_matrix_algorithm(int n, double* A, double* B, double* C, double* D);
void interp_spline(int grid_dim, double step, double* X, double* Y, double* a, double* b, double* c, double* d, double mu1, double mu2);

void get_function(int n, int grid_dim, double& step, double* X, double* Y);
void spline_znachenia(int grid_dim, double* S, double* S2, double* X, double* a, double* b, double* c, double* d);

void spline_pr(int grid_dim, int N, int mu1, int mu2, double* X_control, double* S_control, double* S_control_pr1, double* S_control_pr2, double* a, double* b, double* c, double* d);

void test_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2);
void var1_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2);
void var2_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2);
void var3_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2);
void var4_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2);
void var5_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2);
void var6_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2);
void var7_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2);
void var8_pr(int N, double step2, double* X_control, double* Y_control, double* Y_control_pr1, double* Y_control_pr2);