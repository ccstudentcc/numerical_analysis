#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER

#include "Function.hpp"
#include <stdexcept>
#include <cmath>

class EquationSolver {
protected:
    const Function & F; // Reference to a Function object
public:
    EquationSolver(const Function& F) : F(F) {} // Constructor initializes the function reference
    virtual double solve() = 0; // Pure virtual method for solving the equation
};

// Symbolic function
template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class Bisection_Method : public EquationSolver {
private:
    double a, b; // Interval endpoints
    double eps, delta; // Precision parameters
    int Maxiter; // Maximum iterations
public:
    Bisection_Method(const Function &F, double a, double b, 
        double eps = 1e-9, double delta = 1e-9, int Maxiter = 1000) :
        EquationSolver(F), a(a), b(b), eps(eps), delta(delta), Maxiter(Maxiter) {}
    
    virtual double solve() override {
        // Check if the function has different signs at the endpoints
        if (F(a) * F(b) >= 0) {
            throw std::invalid_argument("Function values at the endpoints must have opposite signs");
        }
        
        double c, u, w, a0, h;
        a0 = a;
        u = sgn(F(a));
        h = b - a;

        for (int iter = 0; iter <= Maxiter; ++iter) {
            h = h / 2.0;
            c = a0 + h; // Calculate midpoint
            if(h < delta || iter == Maxiter-1){
                break;
            }
            w = F(c);

            if (fabs(w) < eps) {
                break; // Found root or reached desired precision
            }
            else if (sgn(w) == u) {
                a0 = c; // Root is in the right half
            }
        }

        return c; // Return the last midpoint
    }
};

class Newton_Method : public EquationSolver {
private:
    double x0; // Initial guess
    double eps; // Precision parameter
    int Maxiter; // Maximum iterations
public:
    Newton_Method(const Function &F, double x0, 
        double eps = 1e-9, int Maxiter = 1000) :
        EquationSolver(F), x0(x0), eps(eps), Maxiter(Maxiter) {}
    
    virtual double solve() override {
        double x = x0;
        for (int iter = 0; iter < Maxiter; ++iter) {
            double fx = F(x); // Function value at x
            if(fabs(fx) < eps){
                break;
            }
            double dfx = F.derivative(x); // Assume F has a derivative method
            if (dfx == 0) {
                throw std::runtime_error("Derivative is zero, cannot proceed");
            }
            
            x = x - fx / dfx; // Newton's method formula
        }
        return x; // Return the last x value
    }
};

class Secant_Method : public EquationSolver {
private:
    double x0, x1; // Initial guess
    double eps, delta; // Precision parameter
    int Maxiter; // Maximum iterations
public:
    Secant_Method(const Function &F, double x0, double x1,
        double eps = 1e-9, double delta = 1e-9, int Maxiter = 1000) :
        EquationSolver(F), x0(x0), x1(x1), eps(eps), delta(delta), Maxiter(Maxiter) {}
    
    virtual double solve() override {
        double xb = x1, xa = x0;
        double u = F(xb), v = F(xa);
        for (int iter = 2; iter < Maxiter; ++iter) {
            double s = (xb - xa) / (u - v);
            xa = xb;
            v = u;
            xb = xb - u * s;

            if(fabs(xb - xa) < delta){
                break;
            }
            u = F(xb);
            if(fabs(u) < eps){
                break;
            }
        }
        return xb; // Return the last x value
    }
};

#endif // EQUATIONSOLVER
