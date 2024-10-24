#include "Function.hpp"
#include "Interpolation.hpp"
#include <iostream>
#include <cmath>
#include <vector>

const double Pi = acos(-1.);

class F1 : public Function {
public:
    double operator() (double x) const {
        return 1.0/(1.0 + 25*pow(x,2));
    }
};


double x_i(int i, int n){
    double x= cos((2.0*i+1.0)/2.0/n*Pi);
    return x;
}

int main() {
    std::cout << "Problem F" << std::endl;
    std::vector<int> n={5,10,15,20};
    F1 f;
    for (auto n_cur : n) {
        std::vector<double> x((size_t)(n_cur), 0.0);
        std::vector<double> y((size_t)(n_cur), 0.0);

        for(int j=0; j<n_cur; j++){ 
            x[j] = x_i(j, n_cur);
            y[j] = f(x[j]);
        }

        NewtonInterpolation func(x,y);
        std::cout << "when n=" << n_cur << ", p(x)=";
        // func.print();
        // std::cout << "Simplify: p(x)=";
        Polynomial func_poly= func.to_polynomial();
        func_poly.print();

    }

    return 0;
}