#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include "Function.hpp"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <numeric> // for std::accumulate
#include <cmath>

/**
 * @brief Represents a point in 2D space.
 */
struct Point
{
    double x, y;

    /**
     * @brief Constructor to initialize a Point.
     *
     * @param x The x-coordinate of the point.
     * @param y The y-coordinate of the point.
     */
    Point(double x = 0, double y = 0) : x(x), y(y) {}

    // Assignment operator
    Point &operator=(const Point &other)
    {
        if (this != &other)
        {
            x = other.x;
            y = other.y;
        }
        return *this;
    }

    // Addition operator
    Point operator+(const Point &other) const
    {
        return Point(x + other.x, y + other.y);
    }

    // Subtraction operator
    Point operator-(const Point &other) const
    {
        return Point(x - other.x, y - other.y);
    }

    // Multiplication operator (element-wise)
    Point operator*(const Point &other) const
    {
        return Point(x * other.x, y * other.y);
    }

    Point operator*(const double &other) const
    {
        return Point(x * other, y * other);
    }

    // Division operator
    Point operator/(const Point &other) const
    {
        if (other.x == 0 || other.y == 0)
        {
            throw std::runtime_error("Division by zero is not allowed.");
        }
        return Point(x / other.x, y / other.y);
    }

    Point operator/(const double &other) const
    {
        if (other == 0)
        {
            throw std::runtime_error("Division by zero is not allowed.");
        }
        return Point(x / other, y / other);
    }

    // Friend function to print the Point
    friend std::ostream &operator<<(std::ostream &os, const Point &point)
    {
        os << point.x << "," << point.y;
        return os;
    }
};

/**
 * @brief Calculate the factorial of n.
 *
 * @param n The number.
 * @return The factorial of n.
 */
int fact(int n){
    if(n==0) return 1;
    return n * fact(n-1);
}

/**
 * @brief Calculate the binomial coefficient "n choose i".
 *
 * @param n The total number of items.
 * @param i The number of items to choose.
 * @return The binomial coefficient.
 */
int binomial(int n, int i)
{
    int res = 1;
    for (int j = 1; j <= i; ++j)
    {
        res *= (n - j + 1) / (double)j;
    }
    return res;
}

/**
 * @brief Abstract base class for interpolation methods.
 */
class Interpolation
{
protected:
    std::vector<double> x_values; ///< Store the x-values of interpolation nodes

public:
    /**
     * @brief Constructor to initialize the x-values.
     *
     * @param x_vals The x-values of interpolation nodes.
     */
    Interpolation(const std::vector<double> &x_vals) : x_values(x_vals) {}

    /**
     * @brief Pure virtual function to evaluate the interpolation polynomial at a given point x.
     *
     * @param x The point at which to evaluate the polynomial.
     * @return The value of the polynomial at x.
     */
    virtual double evaluate(double x) const = 0;

    /**
     * @brief Pure virtual function to print the interpolation polynomial.
     */
    virtual void print() const = 0;
};

/**
 * @brief Class for Newton Interpolation.
 */
class NewtonInterpolation : public Interpolation
{
private:
    std::vector<double> divided_diff; ///< Store divided differences

    /**
     * @brief Helper function to compute divided differences.
     *
     * @param f The function values at the x-values.
     */
    void compute_divided_differences(const std::vector<double> &f)
    {
        int n = x_values.size();
        divided_diff.resize(n);

        // Initialize divided differences with function values
        for (int i = 0; i < n; ++i)
        {
            divided_diff[i] = f[i]; // Assuming f is a vector of function values
        }

        // Compute higher-order divided differences
        for (int j = 1; j < n; ++j)
        {
            for (int i = n - 1; i >= j; --i)
            {
                divided_diff[i] = (divided_diff[i] - divided_diff[i - 1]) / (x_values[i] - x_values[i - j]);
            }
        }
    }

public:
    /**
     * @brief Constructor that takes x-values and function values.
     *
     * @param x_vals The x-values of interpolation nodes.
     * @param f The function values at the x-values.
     */
    NewtonInterpolation(const std::vector<double> &x_vals, const std::vector<double> &f)
        : Interpolation(x_vals)
    {
        compute_divided_differences(f);
    }

    /**
     * @brief Evaluate the Newton interpolation polynomial at a given point x.
     *
     * @param x The point at which to evaluate the polynomial.
     * @return The value of the polynomial at x.
     */
    virtual double evaluate(double x) const override
    {
        int n = x_values.size();
        double result = divided_diff[n - 1]; // Start with the highest divided difference

        // Apply the Newton interpolation formula
        for (int i = n - 2; i >= 0; --i)
        {
            result = result * (x - x_values[i]) + divided_diff[i];
        }

        return result;
    }

    /**
     * @brief Print the Newton interpolation polynomial.
     */
    virtual void print() const override
    {
        int n = x_values.size();

        // Print the first term (constant)
        if (divided_diff[0] != 0)
        {
            std::cout << std::setprecision(7) << divided_diff[0];
        }
        else if (n == 1)
        {
            std::cout << 0 << std::endl;
        }

        // Print the remaining terms
        for (int i = 1; i < n; ++i)
        {
            if (divided_diff[i] != 0.0)
            {
                if (divided_diff[i] > 0.0)
                {
                    if (i == 1 && divided_diff[0] == 0.0)
                    {
                    }
                    else
                    {
                        std::cout << "+";
                    }
                    if (divided_diff[i] != 1.0)
                    {
                        std::cout << std::setprecision(7) << divided_diff[i];
                    }
                }
                else if (divided_diff[i] < 0.0)
                {
                    if (i == 1 && divided_diff[0] == 0.0)
                    {
                    }
                    else
                    {
                        std::cout << "-";
                    }
                    if (divided_diff[i] != -1.0)
                    {
                        std::cout << std::setprecision(7) << -divided_diff[i];
                    }
                }

                // Print the multiplication terms (x - x_0)(x - x_1)...(x - x_{i-1})
                for (int j = 0; j < i; ++j)
                {
                    if (j == 0 && !(divided_diff[i] == 1.0 or divided_diff[i] == -1.0))
                    {
                        std::cout << "*";
                    }
                    else if (0 < j && j< i)
                    {
                        std::cout << "*";
                    }
                    if (x_values[j] == 0.0)
                    {
                        std::cout << "x";
                    }
                    else if (x_values[j] > 0.0)
                    {
                        std::cout << "(x-" << x_values[j] << ")";
                    }
                    else
                    {
                        std::cout << "(x+" << -x_values[j] << ")";
                    }
                }
            }
        }

        std::cout << std::endl;
    }

    /**
     * @brief Convert the Newton interpolation polynomial to a Polynomial object.
     *
     * @return The corresponding Polynomial object.
     */
    Polynomial to_polynomial()
    {
        int n = x_values.size();

        std::vector<double> tem_root({x_values[0]});
        Polynomial f({divided_diff[0]});

        for (int t = 1; t < n; t++)
        {
            Polynomial f_temp({divided_diff[t]});

            for (int j = 0; j < t; j++)
            {
                Polynomial f_multi({-tem_root[j], 1});
                f_temp = f_temp * f_multi;
            }

            f = f + f_temp;

            tem_root.push_back(x_values[t]);
        }

        return f;
    }
};

/**
 * @brief Class for Hermite interpolation.
 *
 * This class implements Hermite interpolation using divided differences
 * to construct the interpolation polynomial.
 */
class HermiteInterpolation : public Interpolation
{
public:
    std::vector<double> y_values; ///< Function values at x values.
    std::vector<int> m_values;    ///< Orders of derivatives at each x value.
    std::vector<double> x_real_values; ///< Real x-values for the divided difference table.
    std::vector<std::vector<double>> divided_diff_table; ///< Divided difference table.
    int total_size; ///< Total size of divided difference table.

    friend std::vector<double> vieteExpansion(const std::vector<double> &roots, const double &an);

    /**
     * @brief Build the divided difference table.
     */
    void build_divided_difference_table()
    {
        total_size = 0;
        for (int mi : m_values)
        {
            total_size += mi + 1; // mi + 1 because we include the function value itself
        }

        divided_diff_table.resize(total_size, std::vector<double>(total_size, 0.0));
        std::vector<std::vector<bool>> flag(total_size, std::vector<bool>(total_size, false));
        std::vector<int> index(x_values.size(), 0);
        x_real_values.resize(total_size, 0.0);

        // Fill the first column with function values and derivatives
        for (size_t i = 0; i < x_values.size(); ++i)
        {
            for (int j = 0; j <= m_values[i]; ++j)
            {
                for (int k = 0; k <= j; ++k)
                {
                    divided_diff_table[index[i] + j][k] = y_values[index[i] + k] / fact(j);
                    flag[index[i] + j][k] = true;
                }
                x_real_values[index[i] + j] = x_values[i];
            }
            if (i + 1 < x_values.size())
            {
                index[i + 1] = index[i] + (m_values[i] + 1);
            }
        }

        // Fill the rest of the table with divided differences
        for (int j = 1; j < total_size; ++j)
        {
            for (int i = j; i < total_size; ++i)
            {
                if (flag[i][j] == false)
                {
                    divided_diff_table[i][j] = (divided_diff_table[i][j - 1] - divided_diff_table[i - 1][j - 1]) /
                                               (x_real_values[i] - x_real_values[i - j]);
                    flag[i][j] = true;
                }
            }
        }
    }

public:
    /**
     * @brief Constructor to initialize Hermite interpolation with data points.
     *
     * @param x_vals The x-values of interpolation nodes.
     * @param y_vals The function values at x-values.
     * @param m_vals The orders of derivatives at each x-value.
     */
    HermiteInterpolation(const std::vector<double> &x_vals, const std::vector<double> &y_vals, const std::vector<int> &m_vals)
        : Interpolation(x_vals), y_values(y_vals), m_values(m_vals)
    {
        build_divided_difference_table();
    }

    /**
     * @brief Evaluate the Hermite polynomial at a given x.
     *
     * @param x The input value.
     * @return The value of the Hermite polynomial at x.
     */
    virtual double evaluate(double x) const override
    {
        int n = x_values.size();
        double result = divided_diff_table[n - 1][n - 1]; // Start with the highest divided difference

        // Apply the Newton interpolation formula
        for (int i = n - 2; i >= 0; --i)
        {
            result = result * (x - x_real_values[i]) + divided_diff_table[i][i];
        }

        return result;
    }

    /**
     * @brief Print the Hermite interpolation polynomial.
     */
    virtual void print() const override
    {
        int n = total_size;

        // Print the first term (constant)
        if (divided_diff_table[0][0] != 0)
        {
            std::cout << std::setprecision(7) << divided_diff_table[0][0];
        }
        else if (n == 1)
        {
            std::cout << 0 << std::endl;
        }

        // Print the remaining terms
        for (int i = 1; i < n; ++i)
        {
            if (divided_diff_table[i][i] != 0.0)
            {
                if (divided_diff_table[i][i] > 0.0)
                {
                    if (i == 1 && divided_diff_table[0][0] == 0.0)
                    {
                    }
                    else
                    {
                        std::cout << "+";
                    }
                    if (divided_diff_table[i][i] != 1.0)
                    {
                        std::cout << std::setprecision(7) << divided_diff_table[i][i];
                    }
                }
                else if (divided_diff_table[i][i] < 0.0)
                {
                    if (i == 1 && divided_diff_table[0][0] == 0.0)
                    {
                    }
                    else
                    {
                        std::cout << "-";
                    }
                    if (divided_diff_table[i][i] != -1.0)
                    {
                        std::cout << std::setprecision(7) << -divided_diff_table[i][i];
                    }
                }

                // Print the multiplication terms (x - x_0)(x - x_1)...(x - x_{i-1})
                for (int j = 0; j < i; ++j)
                {
                    if (j == 0 && !(divided_diff_table[i][i] == 1.0 or divided_diff_table[i][i] == -1.0))
                    {
                        std::cout << "*";
                    }
                    else if (0 < j && j< i)
                    {
                        std::cout << "*";
                    }

                    if (x_real_values[j] == 0.0)
                    {
                        std::cout << "x";
                    }
                    else if (x_real_values[j] > 0.0)
                    {
                        std::cout << "(x-" << x_real_values[j] << ")";
                    }
                    else
                    {
                        std::cout << "(x+" << -x_real_values[j] << ")";
                    }
                }
            }
        }

        std::cout << std::endl;
    }

    /**
     * @brief Convert the Hermite interpolation polynomial to a Polynomial object.
     *
     * @return The corresponding Polynomial object.
     */
    Polynomial to_polynomial()
    {
        int n = total_size;

        std::vector<double> tem_root({x_real_values[0]});
        Polynomial f({divided_diff_table[0][0]});

        for (int t = 1; t < n; t++)
        {
            Polynomial f_temp({divided_diff_table[t][t]});

            for (int j = 0; j < t; j++)
            {
                Polynomial f_multi({-tem_root[j], 1});
                f_temp = f_temp * f_multi;
            }

            f = f + f_temp;

            tem_root.push_back(x_real_values[t]);
        }

        return f;
    }
};

/**
 * @brief Class for Bezier interpolation.
 */
class BezierInterpolation
{
private:
    std::vector<Point> points; ///< Control points for the Bezier curve.
    std::vector<Point> tan_points; ///< Tangent points for each control point.
    const Function &Fx; ///< Function for x-coordinates.
    const Function &Fy; ///< Function for y-coordinates.
    int totalsize; ///< Total number of control points.

public:
    /**
     * @brief Constructor to initialize Bezier interpolation with control points and functions.
     *
     * @param in_points The input points for the Bezier curve.
     * @param Fx The function for x-coordinates.
     * @param Fy The function for y-coordinates.
     */
    BezierInterpolation(const std::vector<Point> &in_points, const Function &Fx, const Function &Fy) 
        : points(in_points), Fx(Fx), Fy(Fy)
    {
        totalsize = points.size();
        tan_points.resize(points.size(), Point(0, 0));

        for (int i = 0; i < totalsize; i++)
        {
            tan_points[i].x = Fx.derivative(points[i].x);
            tan_points[i].y = Fy.derivative(points[i].y);
        }

        for(int i = 0; i < totalsize; i++)
        {
            points[i].x = Fx(points[i].x);
            points[i].y = Fy(points[i].y);
        }
    }

    /**
     * @brief Generate points on the Bezier curve.
     *
     * @param step The step size for interpolation along the curve.
     * @return A vector of points on the Bezier curve.
     */
    std::vector<Point> curve(double step)
    {
        std::vector<Point> regular_point;
        for (int j = 0; j < totalsize - 1; j++)
        {
            std::vector<Point> q((size_t)4, Point(0, 0));
            q[0] = points[j];
            q[1] = points[j];
            q[1].x += tan_points[j].x / sqrt(pow(tan_points[j].x, 2) + pow(tan_points[j].y, 2)) * (points[j + 1].x - points[j].x) / 3;
            q[1].y += tan_points[j].y / sqrt(pow(tan_points[j].x, 2) + pow(tan_points[j].y, 2)) * (points[j + 1].y - points[j].y) / 3;
            q[2] = points[j + 1];
            q[2].x -= tan_points[j].x / sqrt(pow(tan_points[j].x, 2) + pow(tan_points[j].y, 2)) * (points[j + 1].x - points[j].x) / 3;
            q[2].y -= tan_points[j].y / sqrt(pow(tan_points[j].x, 2) + pow(tan_points[j].y, 2)) * (points[j + 1].y - points[j].y) / 3;
            q[3] = points[j + 1];

            int n = 3;
            for (double t = 0; t <= 1; t += step)
            {
                Point res(0, 0);
                for (int i = 0; i <= n; ++i)
                {
                    double b = binomial(n, i) * pow(t, i) * pow(1 - t, n - i);
                    res.x += q[i].x * b;
                    res.y += q[i].y * b;
                }
                regular_point.push_back(res);
            }
        }
        return regular_point;
    }
};

#endif // INTERPOLATION_HPP
