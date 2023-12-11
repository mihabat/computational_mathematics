#include <iostream>
#include <array>
#include <cmath>
#include "lab4.h"
#include <fstream>

double func(double x)
{
    return std::sin(x);
}

int main()
{

    std::ofstream out;
    out.open("5knots.txt");
    out << std::fixed;
    out.precision(16);
    const double I = 1 - std::cos(10);
    for(double h = 0.001; h < 4; h+= 0.00001)
    {
        out << integrate<decltype(func), double, 5>(func, 0, 10, h) - I << std::endl;
    }
    out.close();

    out.open("2knots.txt");
    out << std::fixed;
    out.precision(16);
    for(double h = 0.001; h < 4; h+= 0.00001)
    {
        out << I - integrate<decltype(func), double, 2>(func, 0, 10, h)  << std::endl;
    }
    out.close();
}
