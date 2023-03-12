#include <iostream>
#include <string>
#include <iomanip>
#include "nlohmann/json.hpp"
#include <fstream>
#include <cmath>
#include "Function.h"
//#include <opencv2/opencv.hpp>

using namespace std;
using json = nlohmann::json;

#define PI  3,1415926535 8979323846


int main() {
    std::ifstream l("list.json");
    json data = json::parse(l);

    double step;

    // var //
    // 0 - test
    // 1-4 - main function (Ivlev, Izmaylov, Mitin, Shokurov)
    // 5-8 - function with oscillation (Ivlev, Izmaylov, Mitin, Shokurov) 
    //       Sequence as in the list
    int var;

    int grid_dim; // grid dimension of x
    
    setlocale(LC_ALL,"rus");
    std::cout << "Программа для интерполяции функции с помощью кубических сплайнов" << std::endl;
    std::cout << "Введите размерность сетки: ";
    std::cin >> grid_dim;

    double* X = new double[grid_dim];
    double* Y = new double[grid_dim]; 
    double* Y2 = new double[grid_dim];

    double* a = new double[grid_dim - 1];
    double* b = new double[grid_dim - 1];
    double* c = new double[grid_dim];
    double* d = new double[grid_dim - 1];
    

    std::cout << "Enter variant number: ";
    cin >> var;


    switch (var)
    {
    case 0: {
        std::cout << "---------------" << std::endl;
        std::cout << "Test function !!" << std::endl;
        std::cout << "\tf(x) = " << data["test"]["f_test_1"] << ' ' << ",     x=[-1,0]" << std::endl;
        std::cout << "\tf(x) = " << data["test"]["f_test_2"] << ' ' << ",    x=[0,1]" << std::endl;
        std::cout << "\n\n";

        step = 2.0 / (grid_dim-1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = -1 + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            if (X[i] <= 0)
                Y[i] = X[i] * X[i] * X[i] + 3 * X[i] * X[i];
            else
                Y[i] = -1 * X[i] * X[i] * X[i] + 3 * X[i] * X[i];
        }
        break;
    }
    case 1: {
        std::cout << "---------------" << std::endl;
        std::cout << "Main function !!" << std::endl;
        std::cout << "\tf(x) = " << data["main"]["var_5"]["f_main"] << ' ' << ",     x=[" << data["main"]["var_5"]["a"]<< ',' << data["main"]["var_5"]["b"] << ']' << std::endl;
        std::cout << "\n\n";

        double l = data["main"]["var_5"]["a"];
        double r = data["main"]["var_5"]["b"];

        step = (r - l) / (grid_dim - 1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = log(X[i] + 1) / X[i];
        }

        break;
    }
    case 2:{
        std::cout << "---------------" << std::endl;
        std::cout << "Main function !!" << std::endl;
        std::cout << "\tf(x) = " << data["main"]["var_6"]["f_main"] << ' ' << ",     x=[" << data["main"]["var_6"]["a"] << ',' << data["main"]["var_6"]["b"] << ']' << std::endl;
        std::cout << "\n\n";

        double l = data["main"]["var_6"]["a"];
        double r = data["main"]["var_6"]["b"];

        step = (r-l) / (grid_dim-1);


        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }
        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = sin(X[i] + 1)/X[i];
        }

        for (int i = 0; i < grid_dim; ++i) {
            std::cout <<X[i] << ' ' << Y[i] << std::endl;
        }
        std::cout << "\n\n";
        break;
    }
    case 3: {
        std::cout << "---------------" << std::endl;
        std::cout << "Main function !!" << std::endl;
        std::cout << "\tf(x) = " << data["main"]["var_13"]["f_main"] << ' ' << ",     x=[" << data["main"]["var_13"]["a"] << ',' << data["main"]["var_13"]["b"] << ']' << std::endl;
        std::cout << "\n\n";

        double l = data["main"]["var_13"]["a"];
        double r = data["main"]["var_13"]["b"];

        step = (r - l) / (grid_dim - 1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = (sin(X[i])*sin(X[i])) / X[i];
        }

        break;
    }
    case 4: {
        std::cout << "---------------" << std::endl;
        std::cout << "Main function !!" << std::endl;
        std::cout << "\tf(x) = " << data["main"]["var_20"]["f_main"] << ' ' << ",     x=[" << data["main"]["var_20"]["a"] << ',' << data["main"]["var_20"]["b"] << ']' << std::endl;
        std::cout << "\n\n";

        double l = data["main"]["var_20"]["a"];
        double r = data["main"]["var_20"]["b"];

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
        std::cout << "---------------" << std::endl;
        std::cout << "Oscillating function !!" << std::endl;
        std::cout << "\tf(x) = " << data["main"]["var_5"]["f_oscil"] << ' ' << ",     x=[" << data["main"]["var_5"]["a"] << ',' << data["main"]["var_5"]["b"] << ']' << std::endl;
        std::cout << "\n\n";

        double l = data["main"]["var_5"]["a"];
        double r = data["main"]["var_5"]["b"];

        step = (r - l) / (grid_dim - 1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = log(X[i] + 1) / X[i] + cos(10);
        }

        break;
    }
    case 6: {
        std::cout << "---------------" << std::endl;
        std::cout << "Oscillating function !!" << std::endl;
        std::cout << "\tf(x) = " << data["main"]["var_6"]["f_oscil"] << ' ' << ",     x=[" << data["main"]["var_6"]["a"] << ',' << data["main"]["var_6"]["b"] << ']' << std::endl;
        std::cout << "\n\n";

        double l = data["main"]["var_6"]["a"];
        double r = data["main"]["var_6"]["b"];

        step = (r - l) / (grid_dim - 1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = sin(X[i] + 1) / X[i] + cos(10);
        }

        break;
    }
    case 7: {
        std::cout << "---------------" << std::endl;
        std::cout << "Oscillating function !!" << std::endl;
        std::cout << "\tf(x) = " << data["main"]["var_13"]["f_oscil"] << ' ' << ",     x=[" << data["main"]["var_13"]["a"] << ',' << data["main"]["var_13"]["b"] << ']' << std::endl;
        std::cout << "\n\n";

        double l = data["main"]["var_13"]["a"];
        double r = data["main"]["var_13"]["b"];

        step = (r - l) / (grid_dim - 1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = (sin(X[i])*sin(X[i])) / X[i] + cos(10);
        }

        break;
    }
    case 8: {
        std::cout << "---------------" << std::endl;
        std::cout << "Oscillating function !!" << std::endl;
        std::cout << "\tf(x) = " << data["main"]["var_20"]["f_oscil"] << ' ' << ",     x=[" << data["main"]["var_20"]["a"] << ',' << data["main"]["var_20"]["b"] << ']' << std::endl;
        std::cout << "\n\n";

        double l = data["main"]["var_20"]["a"];
        double r = data["main"]["var_20"]["b"];

        step = (r - l) / (grid_dim - 1);

        for (int i = 0; i < grid_dim; ++i) {
            X[i] = l + i * step;
        }

        for (int i = 0; i < grid_dim; ++i) {
            Y[i] = X[i]*sin(X[i]) /3 + cos(10);
        }

        break;
    }
    }


    interp_spline(grid_dim, step, X, Y, a, b, c, d);

    for (int i = 0; i < grid_dim - 1; i++) {
        std::cout << a[i] << " + ";
        std::cout << b[i] << " * (x-" << X[i + 1] << ") + ";
        std::cout << c[i + 1] / 2 << " * (x-" << X[i + 1] << ")^2 + ";
        std::cout << d[i] / 6 << " * (x-" << X[i + 1] << ")^3" << std::endl;
    }


    std::ofstream Table1("Table_1.csv");

    Table1 << "i, X(i-1), X(i), a(i), b(i), c(i), d(i)\n";
    
    for (int i = 0; i < grid_dim-1; i++) {
        Table1 << i+1 << ',' << X[i + 1] << ',' << X[i] << ',' << a[i] << ',' << b[i] << ',' << c[i+1] << ',' << d[i] << "\n";
    }

    double* S = new double[grid_dim];
    double* S2 = new double[grid_dim];

    S[0] = a[0] + b[0] * (X[0] - X[1]) + c[1] / 2 * (X[0] - X[1]) * (X[0] - X[1]) + d[0] / 6 * (X[0] - X[1]) * (X[0] - X[1]) * (X[0] - X[1]);
    for (int i = 1; i < grid_dim; i++) {
        S[i] = a[i - 1];
    }



    std::ofstream Table2("Table_2.csv");

    Table2 << "j,Xj,F(x),S(x),F(xj)-S(xj),F'(xj),S'(xj),F'(xj)-S'(xj)\n";
    for (int i = 0; i < grid_dim; i++) {
        Table2 << i << ',' << X[i] << ',' << Y[i] << ',' << S[i] << ',' << Y[i]-S[i] << "\n";
    }

    delete[] X;
    delete[] Y;

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

    return 0;
}
