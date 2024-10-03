#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>

const double Pi = acos(-1.);
double l = 89;
double h = 49;
double D = 55;
double beta1 = 11.5 / 180 * Pi;
double A = l*sin(beta1);
double B = l*cos(beta1);
double C = (h+0.5*D)*sin(beta1) - 0.5*D*tan(beta1);
double E = (h+0.5*D)*cos(beta1) - 0.5*D;

class F1 : public Function {
public:
    double operator() (double x) const {
        return A*sin(x)*cos(x) + B*pow(sin(x),2) - C*cos(x) - E*sin(x);
    }
};

void solve_f1() {
    std::cout << "(a) Solved by Newton method," << std::endl;
    Newton_Method solver_f1(F1(), Pi/6);
    double x1 = solver_f1.solve();
    std::cout << "with initial value \\pi/6, the root is: " <<  x1/Pi*180 << "°" << "." <<std::endl;
}

void solve_f2() {
    std::cout << "(b) Solved by Newton method," << std::endl;
    D = 30;
    C = (h+0.5*D)*sin(beta1) - 0.5*D*tan(beta1);
    E = (h+0.5*D)*cos(beta1) - 0.5*D;
    Newton_Method solver_f2(F1(), 33*Pi/180);
    double x1 = solver_f2.solve();
    std::cout << "with initial value 33° but D = 30, the root is: " <<  x1/Pi*180 << "°" << "." <<std::endl;
}

void solve_f3() {
    std::cout << "(c) Solved by secant method." << std::endl;
    Secant_Method solver_f3_1(F1(), 33*Pi/180, 90*Pi/180);
    double x1 = solver_f3_1.solve();
    std::cout << "With initial value 33° and 90°, The root is: " <<  x1/Pi*180 << "°" << "." <<std::endl;

    Secant_Method solver_f3_2(F1(), 33.0*Pi/180.0, 150.0*Pi/180.0);
    double x2 = solver_f3_2.solve();
    std::cout << "With initial value 33° and 150°, The root is: " <<  x2/Pi*180 << "°" << "." <<std::endl;

    Secant_Method solver_f3_3(F1(), 33.0*Pi/180.0, 145.0*Pi/180.0);
    double x3 = solver_f3_3.solve();
    std::cout << "With initial value 33° and 145°, The root is: " <<  x3/Pi*180 << "°" << "." <<std::endl;
}

int main() {
    std::cout << "Problem F" << std::endl;
    solve_f1();
    solve_f2();
    solve_f3();
    return 0;
}