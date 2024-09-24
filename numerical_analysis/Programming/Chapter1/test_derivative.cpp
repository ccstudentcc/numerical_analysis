#include <iostream>
#include "Function.hpp"

int main() {
    // 定义多项式 f(x) = 3x^3 + 2x^2 + x + 5
    Polynomial poly({5, 1, 2, 3});

    double x = 2.0;

    // 计算 f(2)
    std::cout << "f(2) = " << poly(x) << std::endl;

    // 计算 f'(2)
    std::cout << "f'(2) = " << poly.derivative(x) << std::endl;

    // 计算 f''(2)
    std::cout << "f''(2) = " << poly.nth_derivative(x, 2) << std::endl;

    // 计算 f'''(2)
    std::cout << "f'''(2) = " << poly.nth_derivative(x, 3) << std::endl;

    return 0;
}
