#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include "Function.hpp"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>

// Abstract base class for interpolation methods
class Interpolation
{
protected:
    std::vector<double> x_values; // Store the x-values of interpolation nodes
    std::vector<double> coeff;

public:
    // Constructor to initialize the x-values
    Interpolation(const std::vector<double> &x_vals) : x_values(x_vals) {}

    // Pure virtual function to evaluate the interpolation polynomial at a given point x
    virtual double evaluate(double x) const = 0;

    // Pure virtual function to print the interpolation polynomial
    virtual void print() const = 0;
};

// Class for Newton Interpolation
class NewtonInterpolation : public Interpolation
{
private:
    std::vector<double> divided_diff; // Store divided differences

    // Helper function to compute divided differences
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
    // Constructor that takes x-values and function values
    NewtonInterpolation(const std::vector<double> &x_vals, const std::vector<double> &f)
        : Interpolation(x_vals)
    {
        compute_divided_differences(f);
    }

    // Override the evaluate function to compute the Newton interpolation
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

    // Override the print function to print the Newton interpolation polynomial
    virtual void print() const override
    {
        int n = x_values.size();

        // Print the first term (constant)
        if (divided_diff[0] != 0)
        {
            std::cout << std::setprecision(4) << divided_diff[0];
        }

        // Print the remaining terms
        for (int i = 1; i < n; ++i)
        {
            if (divided_diff[i] != 0.0)
            {
                if (divided_diff[i] > 0.0)
                {
                    std::cout << "+";
                    if (divided_diff[i] != 1.0)
                    {
                        std::cout << std::setprecision(4) << divided_diff[i];
                    }
                }
                else if (divided_diff[i] < 0.0)
                {
                    std::cout << "-";
                    if (divided_diff[i] != -1.0)
                    {
                        std::cout << std::setprecision(4) << -divided_diff[i];
                    }
                }

                // Print the multiplication terms (x - x_0)(x - x_1)...(x - x_{i-1})
                for (int j = 0; j < i; ++j)
                {
                    if (j == 0 && !(divided_diff[i] == 1.0 or divided_diff[i] == -1.0))
                    {
                        std::cout << "*";
                    }
                    else if (0 < j < i)
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
     * @brief Convert the Hermite polynomial to its standard form and store the coefficients.
     */
    void print_polynomial()
    {
        // Coefficient vector for the polynomial in normal form
        std::vector<double> coeff(x_values.size(), 0.0);

        // Start with the first divided difference (constant term)
        coeff[0] = divided_diff[0];

        // Build the polynomial using Horner's method
        for (int i = 1; i < x_values.size(); ++i) {
            std::vector<double> temp_coeff(i + 1, 0.0);
            temp_coeff[0] = coeff[0]; // Start with the existing constant term

            // Update temp_coeff by multiplying with (x - x_values[i-1])
            for (int j = 1; j <= i; ++j) {
                temp_coeff[j] = temp_coeff[j - 1] * (-x_values[i - 1]) + (j < coeff.size() ? coeff[j] : 0);
            }

            // Add the divided difference term to the final coefficients
            for (int j = 0; j <= i; ++j) {
                coeff[j] += divided_diff[i] * temp_coeff[j];
            }
        }

        // Create and print the polynomial from coefficients
        Polynomial polynomial(coeff);
        polynomial.print();
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
private:
    std::vector<double> y_values; // Function values at x values. Each f_value is followed by Corresponding derivatives values
    std::vector<int> m_values;    // Orders of derivatives at each x value
    std::vector<double> x_real_values;
    std::vector<std::vector<double>> divided_diff_table; // Divided difference table
    int total_size;                                      // Total size of divided difference table

    // Build the divided difference table
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
                    divided_diff_table[index[i] + j][k] = y_values[index[i] + k];
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
        double result = divided_diff_table[0][n - 1]; // Start with the highest divided difference

        // Apply the Newton interpolation formula
        for (int i = n - 2; i >= 0; --i)
        {
            result = result * (x - x_real_values[i]) + divided_diff_table[i][i];
        }

        return result;
    }

    /**
     * @brief Print the divided difference table.
     */
    virtual void print() const override
    {
        int n = x_values.size();

        // Print the first term (constant)
        if (divided_diff_table[0][0] != 0)
        {
            std::cout << std::setprecision(4) << divided_diff_table[0][0];
        }

        // Print the remaining terms
        for (int i = 1; i < n; ++i)
        {
            if (divided_diff_table[i][i] != 0.0)
            {
                if (divided_diff_table[i][i] > 0.0)
                {
                    std::cout << "+";
                    if (divided_diff_table[i][i] != 1.0)
                    {
                        std::cout << std::setprecision(4) << divided_diff_table[i][i];
                    }
                }
                else if (divided_diff_table[i][i] < 0.0)
                {
                    std::cout << "-";
                    if (divided_diff_table[i][i] != -1.0)
                    {
                        std::cout << std::setprecision(4) << -divided_diff_table[i][i];
                    }
                }

                // Print the multiplication terms (x - x_0)(x - x_1)...(x - x_{i-1})
                for (int j = 0; j < i; ++j)
                {
                    if (j == 0 && !(divided_diff_table[i][i] == 1.0 or divided_diff_table[i][i] == -1.0))
                    {
                        std::cout << "*";
                    }
                    else if (0 < j < i)
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
};

#endif // INTERPOLATION_HPP
