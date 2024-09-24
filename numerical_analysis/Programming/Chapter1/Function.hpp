#ifndef FUNCTION
#define FUNCTION

#include <vector>
#include <cstddef>
#include <stdexcept>

class Function {
public:
    virtual double operator() (double x) const = 0;

    // 返回一阶导数值
    virtual double derivative(double x) const {
        return 0; // 默认实现，如果不重写则返回 0
    }

    // 计算 n 阶导数
    virtual double nth_derivative(double x, int n) const {
        if (n < 0) {
            throw std::invalid_argument("阶数必须是非负整数");
        }
        double result = 0.0;

        // 计算 n 阶导数
        for (int i = 0; i < n; ++i) {
            result = derivative(x); // 每次调用一次导数
            x += 1e-5; // 这里选择一个很小的值来逼近导数
        }
        
        return result;
    }
};

// 示例：实现一个多项式函数 f(x) = a_n * x^n + a_(n-1) * x^(n-1) + ... + a_1 * x + a_0
class Polynomial : public Function {
private:
    std::vector<double> coefficients; // 存储多项式的系数，coefficients[i] 对应 x^i 的系数

public:
    Polynomial(const std::vector<double>& coeffs) : coefficients(coeffs) {}

    // 实现 operator()，计算多项式的值
    virtual double operator() (double x) const override {
        double result = 0.0;
        double power = 1.0; // 从 x^0 开始
        for (size_t i = 0; i < coefficients.size(); ++i) {
            result += coefficients[i] * power;
            power *= x; // 更新 x 的幂次
        }
        return result;
    }

    // 返回当前多项式的导数的多项式
    Polynomial derivative_polynomial() const {
        std::vector<double> deriv_coeffs;
        
        // 如果多项式是常数，导数为 0
        if (coefficients.size() <= 1) {
            deriv_coeffs.push_back(0);
        } else {
            // 计算导数的系数
            for (size_t i = 1; i < coefficients.size(); ++i) {
                deriv_coeffs.push_back(i * coefficients[i]); // a_i * i 对应 x^(i-1)
            }
        }
        
        return Polynomial(deriv_coeffs); // 返回导数多项式
    }

    // 使用 derivative_polynomial 来实现导数计算
    virtual double derivative(double x) const override {
        // 通过导数多项式计算导数值
        Polynomial deriv_poly = this->derivative_polynomial(); // 获取一阶导数多项式
        return deriv_poly(x); // 计算导数多项式在 x 处的值
    }

    // 实现 n 阶导数，返回多项式的 n 阶导数值
    virtual double nth_derivative(double x, int n) const override {
        if (n < 0) {
            throw std::invalid_argument("阶数必须是非负整数");
        }
        
        // 计算第 n 阶导数
        Polynomial current(*this); // 复制当前多项式
        for (int i = 0; i < n; ++i) {
            current = current.derivative_polynomial(); // 取得当前多项式的导数多项式
        }
        return current(x); // 返回第 n 阶导数多项式在 x 处的值
    }
};

#endif
