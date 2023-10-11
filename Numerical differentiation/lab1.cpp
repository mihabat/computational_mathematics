#include <iostream>
#include <array>
#include <cmath>
using namespace std;

template<unsigned int N>
array<double, N + 1> gauss_method(array<array<double, N + 1>, N + 1> A, array<double, N + 1> b)
{
    for(unsigned int i = 0; i < N + 1; i++)      // straight move
    {
        if(A[i][i] == 0)
        {
            for(int j = i + 1; j < N + 1; j++)
            {
                if(A[j][i] != 0)
                {
                    b[i] += b[j];
                    for(int k = 0; k < N + 1; k++)
                    {
                        A[i][k] += A[j][k];
                    }
                    break;
                }
            }
        }
        double a_ii = A[i][i];
        for(int j = i + 1; j < N + 1; j++)
        {
            b[j] = b[j] - b[i] * A[j][i] / a_ii;
            double a_ji = A[j][i];
            for(int k = i; k < N + 1; k++)
            {
                A[j][k] = A[j][k] - A[i][k] * a_ji / a_ii;
            }
        }
        b[i] = b[i] / a_ii;
        for(int j = i; j < N + 1; j++)
        {
            A[i][j] = A[i][j] / a_ii;
        }
    }
    for(int i = N; i > 0; i--)     // reverse move
    {
        for(int j = i - 1; j >= 0; j--)
        {
            b[j] -= b[i] * A[j][i];
            A[j][i] = 0;
        }
    }
    return b;
}

template<typename RealType, unsigned int N>
struct DerivativeCoef
{
    RealType centralCoef;
    array<RealType, N> otherCoefs;
};

template<typename RealType, unsigned int N, unsigned int L>
DerivativeCoef<RealType, N> calcDerivativeCoef(const array<RealType, N>& points) noexcept
{
    array<array<double, N + 1>, N + 1> A;    // A - matrix SLAU
    for(int j = 0; j < N + 1; j++)
        A[0][j] = 1;
    for(int i = 1; i < N + 1; i++)
        A[i][0] = 0;
    int f = 1;
    for(int i = 1; i < N + 1; i++)
    {
        f *= i;  // замена факториалу
        for(int j = 1; j < N + 1; j++)
        {
            A[i][j] = std::pow(points[j - 1], i) / f;
        }
    }
    array<double, N + 1> b;
    for(int i = 0; i < N + 1; i++)
        b[i] = 0;
    b[L] = 1;
    array<double, N + 1> coefs = gauss_method<N>(A, b);
    array<double, N> otherCoef;
    for(int i = 0; i < N; i++)
    {
        otherCoef[i] = coefs[i + 1];
    }
    DerivativeCoef<double, N> S = {coefs[0], otherCoef};
    return S;
}

int main()
{
    const unsigned int N = 5;
    const unsigned int L = 2;
    const array<double, N> points = {-2, -1, 1, 2, 3};
    DerivativeCoef<double, N> S = calcDerivativeCoef<double, N, L>(points);
    cout << fixed;
    cout.precision(16);
    cout << "centerCoef = " << S.centralCoef << " " << endl;
    for(int i = 0; i < N; i++)
    {
        cout << "otherCoef = " << S.otherCoefs[i] << " " << endl;
    }
}
