#include<array>

template<typename xType, typename yType, unsigned int N>
class NewtonInterpolator
{
    private:
        std::array<xType, N> x_values;
        std::array<yType, N> y_values;

    public:
        NewtonInterpolator(const std::array<xType, N> &points, const std::array<yType, N>& values) noexcept
        {
            x_values = points;
            y_values = values;
            for(unsigned int i = 0; i < N; i++)
                for(unsigned int j = N - 1; j > i ; j--)
                    y_values[j] = (y_values[j - 1] - y_values[j]) / (x_values[j - i - 1] - x_values[j]);
        }
        yType interpolate(const xType& x) const noexcept
        {
            yType y_intrpolate = y_values[N - 1];
            for(int i = N - 2; i >= 0; i--)       // GornerScheme
                y_intrpolate = y_values[i] + (x - x_values[i]) * y_intrpolate;

            return y_intrpolate;
        }
};
