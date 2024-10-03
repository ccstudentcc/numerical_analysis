#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>

const double Pi = acos(-1.);

class F1 : public Function {
public:
    double operator() (double x) const {
        return 1.0/x-tan(x);
    }
};

void solve_f1() {
    std::cout << "Solving x^{-1} - \\tan x on [0, \\pi/2]," << std::endl;
    Bisection_Method solver_f1(F1(), 0, Pi/2);
    double x = solver_f1.solve();
    std::cout << "A root is: " << x << "." << std::endl;
}

class F2 : public Function {
public:
    double operator() (double x) const {
        return 1.0/x-pow(2,x);
    }
};

void solve_f2() {
    std::cout << "Solving x^{-1} - 2^x on [0, 1]," << std::endl;
    Bisection_Method solver_f2(F2(), 0, 1);
    double x = solver_f2.solve();
    std::cout << "A root is: " << x << "." << std::endl;
}

class F3 : public Function {
public:
    double operator() (double x) const {
        return pow(2, -x) + exp(x) + 2*cos(x)-6;
    }
};

void solve_f3() {
    std::cout << "Solving 2^{-x} - e^x + 2\\cos x - 6 on [1, 3]," << std::endl;
    Bisection_Method solver_f3(F3(), 1, 3);
    double x = solver_f3.solve();
    std::cout << "A root is: " << x << "." << std::endl;
}

class F4 : public Function {
public:
    double operator() (double x) const {
        return (pow(x, 3) + 4*pow(x, 2) + 3*x + 5)/(2*pow(x, 3) - 9*pow(x, 2) + 18*x - 2);
    }
};

void solve_f4() {
    std::cout << "Solving (x^3+4x^2+3x+5)/(2x^3-9x^2+18-2) on [0, 4]," << std::endl;
    Bisection_Method solver_f4(F4(), 0, 4);
    double x = solver_f4.solve();
    std::cout << "A root is: " << x << "." << std::endl;
}


int main() {
    std::cout << "Problem B" << std::endl;
    solve_f1();
    solve_f2();
    solve_f3();
    solve_f4();
    

    return 0;
}