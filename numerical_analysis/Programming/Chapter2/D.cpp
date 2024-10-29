#include "Function.hpp"
#include "Interpolation.hpp"
#include <iostream>
#include <cmath>
#include <vector>

int main() {
    std::cout << "Problem D" << std::endl;
    std::vector<double> x={0.0,3.0,5.0,8.0,13.0};
    std::vector<double> y={0,75,225,77,383,80,623,74,993,72};
    std::vector<int> m={1,1,1,1,1};

    HermiteInterpolation func(x,y,m);
    std::cout << "(a)";

    Polynomial func_poly(func.to_polynomial());
    
    Polynomial func_deriv = func_poly.derivative_polynomial();
    std::cout << "p(x)=";
    func_poly.print();

    std::cout << "p'(x)=";
    func_deriv.print();

    std::cout << "Speed at t=10s is p'(10)=" << func_deriv(10) << "feet/s. ";
    std::cout << "Position at t=10s is p(10)=" << func_poly(10) << "feet." <<std::endl;
    auto extrema = func_deriv.findExtrema(0.0, 13.0);
    std::cout << "(b) max speed is " << extrema.second.first << "feet/s at t=" << extrema.second.second << "s. ";
    
    std::cout << "it ";
    if(extrema.second.first > 81){
        std::cout <<  "exceeds";
    }
    else{
        std::cout << "doesn't exceed";
    }
    std::cout << " the speed limit." << std::endl;
    



    return 0;
}