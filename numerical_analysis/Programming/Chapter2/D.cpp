#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>

const double Pi = acos(-1.);

class F1 : public Function {
public:
    double operator() (double x) const {
        return sin(x/2) - 1;
    }
};

void solve_f1() {
    std::cout << "Solving \\sin(x/2)-1" << std::endl;
    Secant_Method solver_f1_1(F1(), 0, Pi/2);
    double x1 = solver_f1_1.solve();
    Secant_Method solver_f1_2(F1(), 2*Pi/3, 3*Pi/4);
    double x2 = solver_f1_2.solve();
    std::cout << "With initial values 0 and \\pi/2, the root is: " << x1 << "." << std::endl;
    std::cout << "With initial values 2\\pi/3 and 3\\pi/4, the root is: " << x2 << "." << std::endl;
}

class F2 : public Function {
public:
    double operator() (double x) const {
        return exp(x) - tan(x);
    }
};

void solve_f2() {
    std::cout << "Solving e^x-\\tan(x)" << std::endl;
    Secant_Method solver_f2_1(F2(), 1, 1.4);
    double x1 = solver_f2_1.solve();
    Secant_Method solver_f2_2(F2(), 1.2, 1.6);
    double x2 = solver_f2_2.solve();
    std::cout << "With initial values 1 and 1.4, the root is: " << x1 << "." << std::endl;
    std::cout << "With initial values 1.2 and 1.6, the root is: " << x2 << "." << std::endl;
}

class F3 : public Function {
public:
    double operator() (double x) const {
        return pow(x,3) - 12*pow(x,2) + 3*x + 1;
    }
};

void solve_f3() {
    std::cout << "Solving x^3-12x^2+3x+1" << std::endl;
    Secant_Method solver_f3_1(F3(), 0, -0.5);
    double x1 = solver_f3_1.solve();
    Secant_Method solver_f3_2(F3(), 0.2, -0.7);
    double x2 = solver_f3_2.solve();
    std::cout << "With initial values 0 and -0.5, the root is: " << x1 << "." << std::endl;
    std::cout << "With initial values 0.2 and -0.7, the root is: " << x2 << "." << std::endl;
}

int main() {
    std::cout << "Problem D" << std::endl;
    solve_f1();
    solve_f2();
    solve_f3();

    return 0;
}