#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>

const double Pi = acos(-1.);
const double L = 10;
const double r = 1;
const double V = 12.4; 

class F1 : public Function {
public:
    double operator() (double x) const {
        return L*(0.5*Pi*pow(r,2) - pow(r,2)*asin(x/r) - x*sqrt(pow(r,2)-pow(x,2))) - V;
    }
};

void solve_f1() {
    std::cout << "Solved by Bisection Method," << std::endl;
    Bisection_Method solver_f1(F1(), 0.1, 0.2);
    double x = solver_f1.solve();
    std::cout << "with initial section [0.1,0.2], the depth is: " << x << "ft." << std::endl;
}

void solve_f2() {
    std::cout << "Solved by Newton Method," << std::endl;
    Newton_Method solver_f2(F1(), 0.2);
    double x = solver_f2.solve();
    std::cout << "with initial value 0.2, the depth is: " << x << "ft." << std::endl;
}


void solve_f3() {
    std::cout << "Solved by Secant Method," << std::endl;
    Secant_Method solver_f3(F1(), 0.1, 0.15);
    double x = solver_f3.solve();
    std::cout << "with initial values 0.1 and 0.15 the depth is: " << x << "ft." << std::endl;
}

int main() {
    std::cout << std::fixed << std::setprecision(2);

    std::cout << "Problem E" << std::endl;
    solve_f1();
    solve_f2();
    solve_f3();

    return 0;
}