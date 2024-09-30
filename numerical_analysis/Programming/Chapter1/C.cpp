#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>


class F1 : public Function {
public:
    double operator() (double x) const {
        return x - tan(x);
    }
};

void solve_f1() {
    std::cout << "Solving x = \\tan x near 4.5 and 7.7" << std::endl;
    Newton_Method solver_f1(F1(), 4.5);
    Newton_Method solver_f2(F1(), 7.7);
    double x1 = solver_f1.solve();
    double x2 = solver_f2.solve();
    std::cout << "A root near 4.5 is: " << x1 << std::endl;
    std::cout << "A root near 7.7 is: " << x2 << std::endl;
}

int main() {
    std::cout << "Problem C" << std::endl;
    solve_f1();

    return 0;
}