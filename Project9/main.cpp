#include <iostream>
#include <string>
#include <iomanip>
#include "nlohmann/json.hpp"
#include <fstream>
#include <cmath>
#include "Function.h"
#include "Output.h"

using namespace std;
using json = nlohmann::json;


int main(int argc, char* argv[]) {

    std::ifstream l;
    if (argc > 1)
    {
        l.open(argv[1]);
    }
    else
    {
        l.open("list.json");
    }
    json data = json::parse(l);
    l.close();

    

    int var = data.at("func_id");
    int grid_dim = data.at("n");
    grid_dim++;
    double mu1 = data.at("M1");
    double mu2 = data.at("M2");

    int N = 3 * grid_dim - 2;

    double step;
    double step2;

    double S_Y_error = 0;
    double S_Y_pr1_error = 0;
    double S_Y_pr2_error = 0;

    int X1 = 0;
    int X2 = 0;
    int X3 = 0;

    double X_error1;
    double X_error2;
    double X_error3;
        
    double* a = new double[grid_dim - 1];
    double* b = new double[grid_dim - 1];
    double* c = new double[grid_dim];
    double* d = new double[grid_dim - 1];

    double* X = new double[grid_dim];
    double* Y = new double[grid_dim]; 
    double* S = new double[grid_dim];

    double* X_control = new double[N];
    double* S_control = new double[N];
    double* Y_control = new double[N];

    double* Y_control_pr1 = new double[N];
    double* S_control_pr1 = new double[N];

    double* Y_control_pr2 = new double[N];
    double* S_control_pr2 = new double[N];

    get_function(var, grid_dim, step, X, Y);
    
    step2 = step / 3;

    interp_spline(grid_dim, step, X, Y, a, b, c, d, mu1, mu2);

    switch (var) {
    case 0: {
        test_pr(N, step2, X_control, Y_control, Y_control_pr1, Y_control_pr2);
        break;
    }
    case 1: {
        var1_pr(N, step2, X_control, Y_control, Y_control_pr1, Y_control_pr2);
        break;
    }
    case 2: {
        var2_pr(N, step2, X_control, Y_control, Y_control_pr1, Y_control_pr2);
        break;
    }
    case 3: {
        var3_pr(N, step2, X_control, Y_control, Y_control_pr1, Y_control_pr2);
        break;
    }
    case 4: {
        var4_pr(N, step2, X_control, Y_control, Y_control_pr1, Y_control_pr2);
        break;
    }
    case 5: {
        var5_pr(N, step2, X_control, Y_control, Y_control_pr1, Y_control_pr2);
        break;
    }
    case 6: {
        var6_pr(N, step2, X_control, Y_control, Y_control_pr1, Y_control_pr2);
        break;
    }
    case 7: {
        var7_pr(N, step2, X_control, Y_control, Y_control_pr1, Y_control_pr2);
        break;
    }
    case 8: {
        var8_pr(N, step2, X_control, Y_control, Y_control_pr1, Y_control_pr2);
        break;
    }
    }
    
    
    spline_pr(grid_dim, N, mu1, mu2, X_control, 
        S_control, S_control_pr1, S_control_pr2, a, b, c, d);

    error(N, Y_control, Y_control_pr1, Y_control_pr2, S_control, S_control_pr1, S_control_pr2, S_Y_error, S_Y_pr1_error, S_Y_pr2_error, X1, X2, X3);

    X_error1 = X_control[X1];
    X_error2 = X_control[X2];
    X_error3 = X_control[X3];
    
    std::ofstream out;
    if (argc > 2)
    {
        out.open(argv[2]);
    }
    else
    {
        out.open("Table_1.csv");
    }

    if (!out.is_open())
    {
        throw std::exception();
    }

    Table1(out, grid_dim, X, a, b, c, d);

    out.close();

    if (argc > 3)
    {
        out.open(argv[3]);
    }
    else
    {
        out.open("Table_2.csv");
    }

    if (!out.is_open())
    {
        throw std::exception();
    }


    Table2(out, N, X_control, Y_control, Y_control_pr1, Y_control_pr2, 
        S_control, S_control_pr1, S_control_pr2);



    out.close();

    if (argc > 4) {
        out.open(argv[4]);
    }
    else {
        out.open("Directory.csv");
    }

    if (!out.is_open()) {
        throw std::exception();
    }

    Directory(out, grid_dim, N, S_Y_error, S_Y_pr1_error, S_Y_pr2_error, X_error1, X_error2, X_error3);
    
    out.close();

    //std::cout << 0 << std::endl;
    /*
    for (int i = 0; i < N; i++) {
        double error = Y_control[i] - S_control[i];
        double error2 = Y_control_pr1[i] - S_control_pr1[i];
        double error3 = Y_control_pr2[i] - S_control_pr2[i];
        std::cout << i << ' ' << X_control[i] << ' ' << Y_control[i] << ' ' << S_control[i] << ' ' << error << ' ' << Y_control_pr1[i] << ' ' << S_control_pr1[i] << ' ' << error2 << ' ' << Y_control_pr2[i] << ' ' << S_control_pr2[i] << ' ' << error3 << std::endl;
        //std::cout << Y[0] << std::endl;
        if(abs(error) > S_Y_error)
            S_Y_error = abs(error);
        if (abs(error2 > S_Y_pr1_error))
            S_Y_pr1_error = abs(error2);
        if (abs(error3 > S_Y_pr2_error))
            S_Y_pr2_error = abs(error3);
    }
    */
    //std::cout << "\n\n" << S_Y_error;;

    /*
    for (int i = 0; i < grid_dim - 1; i++) {
        std::cout << a[i] << " + ";
        std::cout << b[i] << " * (x-" << X[i + 1] << ") + ";
        std::cout << c[i + 1] / 2 << " * (x-" << X[i + 1] << ")^2 + ";
        std::cout << d[i] / 6 << " * (x-" << X[i + 1] << ")^3" << std::endl;
    }

    cout << S_Y_error << ' ' << S_Y_pr1_error << ' ' << S_Y_pr2_error << std::endl;
    */
    
    
    

    delete[] X;
    delete[] Y;

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

    delete[] X_control;
    delete[] Y_control;
    delete[] Y_control_pr1;
    delete[] Y_control_pr2;
    delete[] S_control;
    delete[] S_control_pr1;
    delete[] S_control_pr2;

    return 0;
}
