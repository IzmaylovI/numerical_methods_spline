#pragma once
#include <fstream>

void Table1(std::ofstream& out, int grid_dim, double* X, 
	double* a, double* b, double* c, double* d);

void Table2(std::ofstream& out, int N, double* X_control, 
	double* Y_control, double* Y_control_pr1, double* Y_control_pr2, 
	double* S_control, double* S_control_pr1, double* S_control_pr2);

void Directory(std::ofstream& out, double grid_dim, int N, double S_Y_error, double S_Y_pr1_error, double S_Y_pr2_error, double X_error1, double X_error2, double X_error3);