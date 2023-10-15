#include<iostream>
#include<array>
#include "lab2.h"
#define _USE_MATH_DEFINES
#include<cmath>

int main()
{
    std::array<double, 6> b = {2, 1, 0.5, 0.25, 0.125, 0.0625};
    double a = 0;

    std::cout << "N = " << 3 << std::endl;
    for(int i = 0; i < 6; i++)
    {
        const unsigned int N = 3;
        std::array<double, N> noconstX;
        double h = (b[i] - a) / (N - 1);
        for(int j = 0; j < N; j++)
        {
            noconstX[j] = a + j * h;
        }
        const std::array<double, N> points = noconstX;
        std::array<double, N> noconstY;
        for(int j = 0; j < N; j++)
        {
            noconstY[j] = exp(noconstX[j]);
        }
        const std::array<double, N> values = noconstY;
        std::array<double, 1000> unX;
        for(int j = 0; j < 1000; j++)
            unX[j] = a + j * (b[i] - a) / 999;
        const std::array<double, 1000> x = unX;
        double maX = 0;
        NewtonInterpolator<double, double, N> polynom(points, values);
        for(int j = 0; j < 1000; j++)
        {
            if(exp(x[j]) - polynom.interpolate(x[j]) > 0 && exp(x[j]) - polynom.interpolate(x[j]) > maX)
                maX = exp(x[j]) - polynom.interpolate(x[j]);
            if(exp(x[j]) - polynom.interpolate(x[j]) < 0 && polynom.interpolate(x[j])- exp(x[j]) > maX)
                maX = polynom.interpolate(x[j])- exp(x[j]);
        }
        std::cout << std::fixed;
        std::cout.precision(16);
        std::cout << maX << std::endl;
    }

    std::cout << "N = " << 4 << std::endl;
    for(int i = 0; i < 6; i++)
    {
        const unsigned int N = 4;
        std::array<double, N> noconstX;
        double h = (b[i] - a) / (N - 1);
        for(int j = 0; j < N; j++)
        {
            noconstX[j] = a + j * h;
        }
        const std::array<double, N> points = noconstX;
        std::array<double, N> noconstY;
        for(int j = 0; j < N; j++)
        {
            noconstY[j] = exp(noconstX[j]);
        }
        const std::array<double, N> values = noconstY;
        std::array<double, 1000> unX;
        for(int j = 0; j < 1000; j++)
            unX[j] = a + j * (b[i] - a) / 999;
        const std::array<double, 1000> x = unX;
        double maX = 0;
        NewtonInterpolator<double, double, N> polynom(points, values);
        for(int j = 0; j < 1000; j++)
        {
            if(exp(x[j]) - polynom.interpolate(x[j]) > 0 && exp(x[j]) - polynom.interpolate(x[j]) > maX)
                maX = exp(x[j]) - polynom.interpolate(x[j]);
            if(exp(x[j]) - polynom.interpolate(x[j]) < 0 && polynom.interpolate(x[j])- exp(x[j]) > maX)
                maX = polynom.interpolate(x[j])- exp(x[j]);
        }
        std::cout << std::fixed;
        std::cout.precision(16);
        std::cout << maX << std::endl;
    }

    std::cout << "N = " << 5 << std::endl;
    for(int i = 0; i < 6; i++)
    {
        const unsigned int N = 5;
        std::array<double, N> noconstX;
        double h = (b[i] - a) / (N - 1);
        for(int j = 0; j < N; j++)
        {
            noconstX[j] = a + j * h;
        }
        const std::array<double, N> points = noconstX;
        std::array<double, N> noconstY;
        for(int j = 0; j < N; j++)
        {
            noconstY[j] = exp(noconstX[j]);
        }
        const std::array<double, N> values = noconstY;
        std::array<double, 1000> unX;
        for(int j = 0; j < 1000; j++)
            unX[j] = a + j * (b[i] - a) / 999;
        const std::array<double, 1000> x = unX;
        double maX = 0;
        NewtonInterpolator<double, double, N> polynom(points, values);
        for(int j = 0; j < 1000; j++)
        {
            if(exp(x[j]) - polynom.interpolate(x[j]) > 0 && exp(x[j]) - polynom.interpolate(x[j]) > maX)
                maX = exp(x[j]) - polynom.interpolate(x[j]);
            if(exp(x[j]) - polynom.interpolate(x[j]) < 0 && polynom.interpolate(x[j])- exp(x[j]) > maX)
                maX = polynom.interpolate(x[j])- exp(x[j]);
        }
        std::cout << std::fixed;
        std::cout.precision(16);
        std::cout << maX << std::endl;
    }

    std::cout << "N = " << 3 << "(Chebyshev's knots)" << std::endl;
    for(int i = 0; i < 6; i++)
    {
        const unsigned int N = 3;
        std::array<double, N> noconstX;
        for(int j = 0; j < N; j++)
        {
            noconstX[j] = (b[i] + a) / 2 + (b[i] - a) / 2 * cos(M_PI * (2 * j + 1) / (2 * N));
        }
        const std::array<double, N> points = noconstX;
        std::array<double, N> noconstY;
        for(int j = 0; j < N; j++)
        {
            noconstY[j] = exp(noconstX[j]);
        }
        const std::array<double, N> values = noconstY;
        std::array<double, 1000> unX;
        for(int j = 0; j < 1000; j++)
            unX[j] = a + j * (b[i] - a) / 999;
        const std::array<double, 1000> x = unX;
        double maX = 0;
        NewtonInterpolator<double, double, N> polynom(points, values);
        for(int j = 0; j < 1000; j++)
        {
            if(exp(x[j]) - polynom.interpolate(x[j]) > 0 && exp(x[j]) - polynom.interpolate(x[j]) > maX)
                maX = exp(x[j]) - polynom.interpolate(x[j]);
            if(exp(x[j]) - polynom.interpolate(x[j]) < 0 && polynom.interpolate(x[j])- exp(x[j]) > maX)
                maX = polynom.interpolate(x[j])- exp(x[j]);
        }
        std::cout << std::fixed;
        std::cout.precision(16);
        std::cout << maX << std::endl;
    }

    std::cout << "N = " << 4 << "(Chebyshev's knots)" << std::endl;
    for(int i = 0; i < 6; i++)
    {
        const unsigned int N = 4;
        std::array<double, N> noconstX;
        for(int j = 0; j < N; j++)
        {
            noconstX[j] = (b[i] + a) / 2 + (b[i] - a) / 2 * cos(M_PI * (2 * j + 1) / (2 * N));
        }
        const std::array<double, N> points = noconstX;
        std::array<double, N> noconstY;
        for(int j = 0; j < N; j++)
        {
            noconstY[j] = exp(noconstX[j]);
        }
        const std::array<double, N> values = noconstY;
        std::array<double, 1000> unX;
        for(int j = 0; j < 1000; j++)
            unX[j] = a + j * (b[i] - a) / 999;
        const std::array<double, 1000> x = unX;
        double maX = 0;
        NewtonInterpolator<double, double, N> polynom(points, values);
        for(int j = 0; j < 1000; j++)
        {
            if(exp(x[j]) - polynom.interpolate(x[j]) > 0 && exp(x[j]) - polynom.interpolate(x[j]) > maX)
                maX = exp(x[j]) - polynom.interpolate(x[j]);
            if(exp(x[j]) - polynom.interpolate(x[j]) < 0 && polynom.interpolate(x[j])- exp(x[j]) > maX)
                maX = polynom.interpolate(x[j])- exp(x[j]);
        }
        std::cout << std::fixed;
        std::cout.precision(16);
        std::cout << maX << std::endl;
    }

    std::cout << "N = " << 5 << "(Chebyshev's knots)" << std::endl;
    for(int i = 0; i < 6; i++)
    {
        const unsigned int N = 5;
        std::array<double, N> noconstX;
        for(int j = 0; j < N; j++)
        {
            noconstX[j] = (b[i] + a) / 2 + (b[i] - a) / 2 * cos(M_PI * (2 * j + 1) / (2 * N));
        }
        const std::array<double, N> points = noconstX;
        std::array<double, N> noconstY;
        for(int j = 0; j < N; j++)
        {
            noconstY[j] = exp(noconstX[j]);
        }
        const std::array<double, N> values = noconstY;
        std::array<double, 1000> unX;
        for(int j = 0; j < 1000; j++)
            unX[j] = a + j * (b[i] - a) / 999;
        const std::array<double, 1000> x = unX;
        double maX = 0;
        NewtonInterpolator<double, double, N> polynom(points, values);
        for(int j = 0; j < 1000; j++)
        {
            if(exp(x[j]) >= polynom.interpolate(x[j]) && exp(x[j]) - polynom.interpolate(x[j]) > maX)
                maX = exp(x[j]) - polynom.interpolate(x[j]);
            if(exp(x[j]) < polynom.interpolate(x[j]) && polynom.interpolate(x[j])- exp(x[j]) > maX)
                maX = polynom.interpolate(x[j])- exp(x[j]);
        }
        std::cout << std::fixed;
        std::cout.precision(16);
        std::cout << maX << std::endl;
    }
}
