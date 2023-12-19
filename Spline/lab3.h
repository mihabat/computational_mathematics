#include <type_traits>
#include <vector>
#include <array>

template <typename RealType>
int find_lower_boundary(const std::vector<RealType> &array, RealType x) {
    int low = 0;
    int high = array.size() - 1;
    int mid;
    while (low <= high) {
        mid = low + (high - low) / 2;
        if (array[mid] == x) {
            return mid;
        } else if (array[mid] < x) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return low - 1;
}

template <typename Type> class ThreeDiagonalMatrix {
private:
    std::vector<std::array<Type, 3>> m_data;

public:
    ThreeDiagonalMatrix(const std::vector<std::array<Type, 3>> &data)
        : m_data{data} {}
    Type operator()(const unsigned int i, const unsigned int j) const {
        if (j == i - 1) {
            return m_data[i][0];
        } else if (j == i) {
            return m_data[i][1];
        } else if (j == i + 1) {
            return m_data[i][2];
        } else {
            return 0;
        }
    }
};

template <typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

template <typename Type>
using DiffType = decltype(std::declval<Type>() - std::declval<Type>());

template <typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve(const ThreeDiagonalMatrix<mType> &matrix, const std::vector<cType> &column)
{
    std::vector<double> p_vector{0}, q_vector{0};
    unsigned int N = column.size() - 1;
    p_vector.reserve(N);
    q_vector.reserve(N);
    for (int i = 0; i < N + 1; ++i)
    {
        p_vector.push_back(-1 * matrix(i, i + 1) / (matrix(i, i - 1) * p_vector[i] + matrix(i, i)));
        q_vector.push_back((column[i] - matrix(i, i - 1) * q_vector[i]) / (matrix(i, i - 1) * p_vector[i] + matrix(i, i)));
    }
    std::vector<double> solution(N + 1);
    solution[N] = (column[N] - matrix(N, N - 1) * q_vector[N]) / (matrix(N, N - 1) * p_vector[N] + matrix(N, N));
    for (int i = N - 1; i >= 0; --i)
        solution[i] = p_vector[i + 1] * solution[i + 1] + q_vector[i + 1];
    return solution;
}
template <typename xType, typename yType> class CubicSpline
{
private:
    std::vector<xType> m_points;
    std::vector<yType> m_values;
    std::vector<DivisType<yType, xType>> m_a;
    std::vector<DivisType<yType, xType>> m_b;

public:
    CubicSpline(const std::vector<xType> &points, const std::vector<yType> &values): m_points{points}, m_values{values}
    {
        unsigned int N = points.size();
        std::vector<xType> points_differences;
        points_differences.reserve(N - 1);
        for (int i = 0; i < N - 1; ++i)
            points_differences.push_back(points[i + 1] - points[i]);
        std::vector<yType> values_differences;
        values_differences.reserve(N - 1);
        for (int i = 0; i < N - 1; ++i)
            values_differences.push_back(values[i + 1] - values[i]);

        std::vector<std::array<DivisType<int, xType>, 3>> data;
        data.reserve(N);
        data.push_back({0, 2 / points_differences[0], 1 / points_differences[0]});
        for (unsigned int i = 0; i < N - 2; ++i)
        {
            data.push_back(
                {1 / points_differences[i],
                 (2 / points_differences[i] + 2 / points_differences[i + 1]),
                 1 / points_differences[i + 1]});
        }
        data.push_back({1 / points_differences[N - 2], 2 / points_differences[N - 2], 0});
        ThreeDiagonalMatrix matrix(data);

        std::vector<DivisType<yType, xType>> column;
        column.reserve(N);
        column.push_back(3 * (values_differences[0]) /
                         (points_differences[0] * points_differences[0]));
        for (int i = 0; i < N - 2; ++i)
        {
            column.push_back(
                3 * (values_differences[i] /
                         (points_differences[i] * points_differences[i]) +
                     values_differences[i + 1] / (points_differences[i + 1] *
                                                  points_differences[i + 1])));
        }
        column.push_back(
            3 * (values_differences[N - 2]) /
            (points_differences[N - 2] * points_differences[N - 2]));

        std::vector<DivisType<yType, xType>> k = solve(matrix, column);

        m_a.reserve(N - 1);
        for (unsigned int i = 0; i < N - 1; ++i)
        {
            m_a.push_back(k[i] * points_differences[i] - (values_differences[i]));
        }

        m_b.reserve(N - 1);
        for (unsigned int i = 0; i < N - 1; ++i)
        {
            m_b.push_back((-1) * k[i + 1] * points_differences[i] + (values_differences[i]));
        }
    }

    yType interpolate(const xType &x) const noexcept
    {
        int k = find_lower_boundary(m_points, x);
        xType t = (x - m_points[k]) / (m_points[k + 1] - m_points[k]);
        return (1 - t) * m_values[k] + t * m_values[k + 1] + t * (1 - t) * ((1 - t) * m_a[k] + t * m_b[k]);
    };
};
