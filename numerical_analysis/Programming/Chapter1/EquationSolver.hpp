#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER

#include "Function.hpp"
#include <stdexcept>
#include <cmath>

class EquationSolver {
protected:
    const Function & F;
public:
    EquationSolver(const Function& F) : F(F) {}
    virtual double solve() = 0;
};

class Bisection_Method : public EquationSolver {
private:
    double a, b;
    double eps, delta;
    int Maxiter;
public:
    Bisection_Method(const Function &F, double a, double b, 
        double eps = 1e-7, double delta = 1e-6, int Maxiter = 50) :
        EquationSolver(F), a(a), b(b), eps(eps), delta(delta), Maxiter(Maxiter) {}
    
    virtual double solve() override {
        if (F(a) * F(b) >= 0) {
            throw std::invalid_argument("函数在区间端点的值必须有不同符号");
        }
        
        double c;
        for (int iter = 0; iter < Maxiter; ++iter) {
            c = (a + b) / 2.0; // 中点
            if (F(c) == 0.0 || (b - a) / 2.0 < eps) {
                return c; // 找到根或达到精度
            }
            
            if (F(c) * F(a) < 0) {
                b = c; // 根在左侧
            } else {
                a = c; // 根在右侧
            }
        }
        return c; // 返回最后的中点
    }
};

class Newton_Method : public EquationSolver {
private:
    double x0;
    double eps;
    int Maxiter;
public:
    Newton_Method(const Function &F, double x0, 
        double eps = 1e-7, int Maxiter = 50) :
        EquationSolver(F), x0(x0), Maxiter(Maxiter), eps(eps) {}
    
    virtual double solve() override {
        double x = x0;
        for (int iter = 0; iter < Maxiter; ++iter) {
            double fx = F(x);
            double dfx = F.derivative(x); // 假设 F 有导数方法
            if (dfx == 0) {
                throw std::runtime_error("导数为零，无法继续");
            }
            
            double x_new = x - fx / dfx;
            if (std::abs(x_new - x) < eps) {
                return x_new; // 找到根或达到精度
            }
            x = x_new; // 更新 x
        }
        return x; // 返回最后的 x 值
    }
};

#endif
