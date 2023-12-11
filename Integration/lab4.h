#include <array>
#include <type_traits>

template<typename A>
struct ArgumentGetter;

template <typename RealType, std::size_t N> struct Knots;
template <typename RealType> struct Knots<RealType, 5>
{
    static constexpr std::array<RealType, 5> knots{-0.9061798459386640, -0.5384693101056831, 0, 0.5384693101056831, 0.9061798459386640};
};
template <typename RealType> struct Knots<RealType, 2>
{
    static constexpr std::array<RealType, 2> knots{-0.5773502691896258, 0.5773502691896258};
};

template <typename RealType, std::size_t N> struct Weights;
template <typename RealType> struct Weights<RealType, 5>
{
    static constexpr std::array<RealType, 5> weights{0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891};
};
template <typename RealType> struct Weights<RealType, 2>
{
    static constexpr std::array<RealType, 2> weights{1, 1};
};

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename T>
using Dif = decltype(std::declval<T>() - std::declval<T>());

/* Функция производит интегрирование на одном отрезке */
template<typename Callable, typename RealType, std::size_t N>
decltype(auto) integrate(
    const Callable& func,  // Интегрируемая функция
    const typename ArgumentGetter<Callable>::Argument& start,  // начало отрезка
    const typename ArgumentGetter<Callable>::Argument& end)  // конец отрезка
{
    RealType I = 0;
    for (int i = 0; i < N; i++)
    {
        I += Weights<RealType, N>::weights[i] * func((end + start) / 2 + (end - start) / 2 * Knots<RealType, N>::knots[i]);
    }
    return (end - start) / 2 * I;
}

/* Функция производит интегрирование, разбивая отрезок на подотрезки длиной не более dx */
template<typename Callable, typename RealType, std::size_t N>
decltype(auto) integrate(
    const Callable& func,  // Интегрируемая функция
    const typename ArgumentGetter<Callable>::Argument& start,  // начало отрезка
    const typename ArgumentGetter<Callable>::Argument& end,  // конец отрезка
    const Dif<typename ArgumentGetter<Callable>::Argument>& dx)  // Длина подотрезка
{
    unsigned int Count = (end - start) / dx + 1;
    double h = (end - start) / Count;
    RealType I = 0;
    for (int i = 0; i < Count; i++)
    {
        I += integrate<Callable, RealType, N>(func, start + i * h, start + (i + 1) * h);
    }
    return I;
}
