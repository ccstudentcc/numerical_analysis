#ifndef FUNCTION
#define FUNCTION

#include <vector>
#include <cstddef>
#include <stdexcept>

class Function {
public:
    virtual double operator() (double x) const = 0;

    // Return the value of the first derivative using central difference
    virtual double derivative(double x) const {
        const double h = 1e-8; // A small increment
        return ((*this)(x + h) - (*this)(x - h)) / (2 * h); // Central difference formula
    }

    // Calculate the n-th derivative
    virtual double nth_derivative(double x, int n) const {
        return 0; // Default implementation returns 0 if not overridden
    }
};

// Example: Implement a polynomial function f(x) = a_n * x^n + a_(n-1) * x^(n-1) + ... + a_1 * x + a_0
class Polynomial : public Function {
private:
    std::vector<double> coefficients; // Store the coefficients of the polynomial, coefficients[i] corresponds to the coefficient of x^i

public:
    Polynomial(const std::vector<double>& coeffs) : coefficients(coeffs) {}

    // Implement operator() to calculate the value of the polynomial
    virtual double operator() (double x) const override {
        double result = 0.0;
        double power = 1.0; // Start from x^0
        for (size_t i = 0; i < coefficients.size(); ++i) {
            result += coefficients[i] * power;
            power *= x; // Update the power of x
        }
        return result;
    }

    // Return the derivative as a polynomial
    Polynomial derivative_polynomial() const {
        std::vector<double> deriv_coeffs;
        
        // If the polynomial is a constant, the derivative is 0
        if (coefficients.size() <= 1) {
            deriv_coeffs.push_back(0);
        } else {
            // Calculate the coefficients of the derivative
            for (size_t i = 1; i < coefficients.size(); ++i) {
                deriv_coeffs.push_back(i * coefficients[i]); // a_i * i corresponds to x^(i-1)
            }
        }
        
        return Polynomial(deriv_coeffs); // Return the derivative polynomial
    }

    // Use derivative_polynomial to implement derivative calculation
    virtual double derivative(double x) const override {
        // Calculate the derivative value using the derivative polynomial
        Polynomial deriv_poly = this->derivative_polynomial(); // Get the first derivative polynomial
        return deriv_poly(x); // Calculate the value of the derivative polynomial at x
    }

    // Implement n-th derivative, returning the value of the polynomial's n-th derivative
    virtual double nth_derivative(double x, int n) const override {
        if (n < 0) {
            throw std::invalid_argument("Order must be a non-negative integer");
        }
        
        // Calculate the n-th derivative
        Polynomial current(*this); // Copy the current polynomial
        for (int i = 0; i < n; ++i) {
            current = current.derivative_polynomial(); // Get the derivative polynomial of the current polynomial
        }
        return current(x); // Return the value of the n-th derivative polynomial at x
    }
};

#endif // FUNCTION
