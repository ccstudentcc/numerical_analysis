#ifndef FUNCTION
#define FUNCTION

#include <vector>
#include <cstddef>
#include <stdexcept>
#include <iostream>
#include <iomanip> // For std::setprecision

/**
 * @brief Abstract base class for mathematical functions.
 * 
 * This class provides a virtual interface for functions that can be evaluated,
 * differentiated, and have their n-th derivatives calculated.
 */
class Function {
public:
    /**
     * @brief Evaluate the function at a given x.
     * 
     * @param x The input value.
     * @return The value of the function at x.
     */
    virtual double operator() (double x) const = 0;

    /**
     * @brief Calculate the first derivative of the function at a given x using central difference.
     * 
     * @param x The point at which to evaluate the derivative.
     * @return The value of the first derivative at x.
     */
    virtual double derivative(double x) const {
        const double h = 1e-8; // A small increment
        return ((*this)(x + h) - (*this)(x - h)) / (2 * h); // Central difference formula
    }

    /**
     * @brief Calculate the n-th derivative of the function at a given x.
     * 
     * @param x The point at which to evaluate the n-th derivative.
     * @param n The order of the derivative.
     * @return The value of the n-th derivative at x.
     * @throws std::invalid_argument if n is negative.
     */
    virtual double nth_derivative(double x, int n) const {
        return 0; // Default implementation returns 0 if not overridden
    }
};

/**
 * @brief Class for representing polynomial functions.
 * 
 * This class implements a polynomial function of the form
 * f(x) = a_n * x^n + a_(n-1) * x^(n-1) + ... + a_1 * x + a_0.
 */
class Polynomial : public Function {
private:
    std::vector<double> coefficients; ///< Store the coefficients of the polynomial

public:
    /**
     * @brief Constructor that initializes the polynomial coefficients.
     * 
     * @param coeffs A vector of coefficients for the polynomial. Defaults to {0.0}.
     */
    Polynomial(const std::vector<double>& coeffs={0.0}) : coefficients(coeffs) {}

    /**
     * @brief Copy constructor.
     * 
     * @param other The Polynomial object to copy.
     */
    Polynomial(const Polynomial& other) : coefficients(other.coefficients) {}

    /**
     * @brief Copy assignment operator.
     * 
     * @param other The Polynomial object to copy.
     * @return A reference to the current object.
     */
    Polynomial& operator=(const Polynomial& other) {
        if (this != &other) { // Check for self-assignment
            coefficients = other.coefficients; // Copy coefficients
        }
        return *this; // Return the current object
    }

    /**
     * @brief Evaluate the polynomial at a given x.
     * 
     * @param x The input value.
     * @return The value of the polynomial at x.
     */
    virtual double operator() (double x) const override {
        double result = 0.0;
        double power = 1.0; // Start from x^0
        for (size_t i = 0; i < coefficients.size(); ++i) {
            result += coefficients[i] * power;
            power *= x; // Update the power of x
        }
        return result;
    }

    /**
     * @brief Print the polynomial in human-readable form.
     */
    void print() const {
        int n=coefficients.size();

        // Print the first term (constant)
        if(coefficients[0]!=0.0){
            std::cout << std::setprecision(4) << coefficients[0];
        }
        

        // Print the remaining terms
        for (int i = 1; i < n; ++i) {
            if (coefficients[i] != 0.0) {
                if (coefficients[i]>0.0) {
                    std::cout << "+";
                    if(coefficients[i]!=1.0){
                        std::cout << std::setprecision(4) << coefficients[i];
                    }
                }
                else if(coefficients[i]<0.0){
                    std::cout << "-";
                    if(coefficients[i]!=-1.0){
                        std::cout << std::setprecision(4) << -coefficients[i];
                    }
                }
                if(!(coefficients[i]==1.0 or coefficients[i]==-1.0)){
                    std::cout << "*";
                }

                std::cout << "x^" << i;

            }
        }
        std::cout << std::endl;
    }

    /**
     * @brief Return the first derivative of the polynomial as a new Polynomial.
     * 
     * @return A Polynomial object representing the derivative.
     */
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

    /**
     * @brief Calculate the value of the first derivative at a given x.
     * 
     * @param x The point at which to evaluate the derivative.
     * @return The value of the derivative at x.
     */
    virtual double derivative(double x) const override {
        // Calculate the derivative value using the derivative polynomial
        Polynomial deriv_poly = this->derivative_polynomial(); // Get the first derivative polynomial
        return deriv_poly(x); // Calculate the value of the derivative polynomial at x
    }

    /**
     * @brief Calculate the n-th derivative of the polynomial and return it as a new Polynomial.
     * 
     * This function computes the n-th derivative of the current polynomial
     * and returns it as a new Polynomial object.
     * 
     * @param n The order of the derivative to calculate. Must be a non-negative integer.
     * @return A Polynomial object representing the n-th derivative.
     * @throws std::invalid_argument if n is negative.
     */
    Polynomial nth_derivative_polynomial(int n) const {
        if (n < 0) {
            throw std::invalid_argument("Order must be a non-negative integer");
        }

        Polynomial current(*this); ///< Create a copy of the current polynomial
        for (int i = 0; i < n; ++i) {
            current = current.derivative_polynomial(); ///< Get the derivative polynomial of the current polynomial
        }
        return current; ///< Return the n-th derivative polynomial
    }

    /**
     * @brief Calculate the value of the n-th derivative at a given point.
     * 
     * This function calculates the n-th derivative of the polynomial
     * and evaluates it at the specified point x.
     * 
     * @param x The point at which to evaluate the n-th derivative.
     * @param n The order of the derivative to calculate. Must be a non-negative integer.
     * @return The value of the n-th derivative at point x.
     * @throws std::invalid_argument if n is negative.
     */
    virtual double nth_derivative(double x, int n) const override {
        if (n < 0) {
            throw std::invalid_argument("Order must be a non-negative integer");
        }

        Polynomial nth_deriv_poly = this->nth_derivative_polynomial(n); ///< Get the n-th derivative polynomial
        return nth_deriv_poly(x); ///< Return the value of the n-th derivative polynomial at x
    }
};


#endif // FUNCTION
