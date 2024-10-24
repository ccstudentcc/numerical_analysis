#include "Function.hpp"
#include "Interpolation.hpp"
#include <iostream>
#include <cmath>
#include <vector>

int main() {
    std::cout << "Problem E" << std::endl;
    std::vector<double> day={0,6,10,13,17,20,28};
    std::vector<double> sp1={6.67,17.3,42.7,37.3,30.1,29.3,28.7};
    std::vector<double> sp2={6.67,16.1,18.9,15.0,10.6,9.44,8.89};

        NewtonInterpolation func1(day, sp1);
        Polynomial func_poly1= func1.to_polynomial();
        std::cout << "(a)For sp1" << ", sp1(x)=";
        func_poly1.print();

        NewtonInterpolation func2(day, sp2);
        Polynomial func_poly2= func2.to_polynomial();
        std::cout << "For sp2" << ", sp2(x)=";
        func_poly2.print();

        std::cout << "(b)For sp1" << ", sp1(43)=" << func_poly1(43) << std::endl;

        std::cout << "For sp2" << ", sp2(43)=" << func_poly2(43) << std::endl;

        

    return 0;
}