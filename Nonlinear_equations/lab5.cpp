#define _USE_MATH_DEFINES
#define M_PI

#include "lab5.h"
#include <iostream>
#include <array>
#include <fstream>
#include <string>

double func(double x)
{
    return (x * x + std::tan(x) * std::tan(x) - 1);
}


int main()
{
    std::array<double, 4> ecc = {0.1, 0.2, 0.5, 0.8};
    std::array<std::string, 4> name = {"e1.txt", "e2.txt", "e3.txt", "e4.txt"};
    std::ofstream out;
    out << std::fixed;
    out.precision(16);
    double tol = 0.00001;
    for(unsigned int i = 0; i < 4; i++)
    {
        out.open(name[i]);
        for(unsigned int j = 0; j < 50; j++)
        {
            out << keplerSolver(ecc[i], M_PI / 4., j, tol) << std::endl;
        }
        out.close();
    }

    double t_opt1 = -2 / (2 * 0.7 + 2 * std::tan(0.7) / (std::cos(0.7) * std::cos(0.7)) + 2 * 0.6 + 2 * std::tan(0.6) / (std::cos(0.6) * std::cos(0.6)));
    double t_opt2 = 2 / (2 * 0.7 + 2 * std::tan(0.7) / (std::cos(0.7) * std::cos(0.7)) + 2 * 0.6 + 2 * std::tan(0.6) / (std::cos(0.6) * std::cos(0.6)));
    const double tau1 = t_opt1;  // t = t_opt
    const double tau2 = t_opt2;  // t = t_opt
    double x1 = solve(func, tau1, 0.7, 20);
    double x2 = solve(func, tau2, -0.7, 20);
    double y1 = std::tan(x1);
    double y2 = std::tan(x2);

    std::cout << std::fixed;
    std::cout.precision(6);
    std::cout << "x1 = " << x1 << " y1 = " << y1 << std::endl;
    std::cout << "x2 = " << x2 << " y2 = " << y2;
}
