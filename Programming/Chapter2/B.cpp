#include "Function.hpp"
#include "Interpolation.hpp"
#include <iostream>
#include <cmath>
#include <vector>

class F1 : public Function {
public:
    double operator() (double x) const {
        return 1.0/(1.0 + pow(x,2));
    }
};


double x_i(int i, int n){
    double x= -5.0 + 10.0*i/n;
    return x;
}

int main() {
    std::cout << "Problem B" << std::endl;
    std::vector<int> n={2,4,6,8};
    F1 f;
    for (auto n_cur : n) {
        std::vector<double> x((size_t)(n_cur+1), 0.0);
        std::vector<double> y((size_t)(n_cur+1), 0.0);

        for(int j=0; j<=n_cur; j++){ 
            x[j] = x_i(j, n_cur);
            y[j] = f(x[j]);
        }

        NewtonInterpolation func(x,y);
        Polynomial func_poly= func.to_polynomial();
        std::cout << "when n=" << n_cur << ", p(x)=";
        func_poly.print();
    }

    return 0;
}