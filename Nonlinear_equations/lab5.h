#include <cmath>

double keplerSolver(double ecc, double meanAnomaly, unsigned int maxIter, double tol) // Newton's method
{
    double E = meanAnomaly; //
    double currentE = E;
    for(unsigned int i = 0; i < maxIter; i++)
    {
        currentE = E - (E - ecc * std::sin(E) - meanAnomaly) / (1 - ecc * std::cos(E));
        if(std::abs(currentE - E) <= tol)
            break;
        E = currentE;
    }
    return currentE;
}


template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename Callable, typename RealType>
decltype(auto) solve(
    const Callable& func,
    const RealType& tau,
    const typename ArgumentGetter<Callable>::Argument& initialGuess,
    const unsigned int nIteration
                    )
{
    typename ArgumentGetter<Callable>::Argument E = initialGuess;
    for(unsigned int i = 0; i < nIteration; i++)
        E = E + tau * func(E);
    return E;
}
