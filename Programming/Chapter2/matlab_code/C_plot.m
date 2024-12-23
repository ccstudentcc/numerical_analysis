% C_plot.m

clear;

% 符号工具箱
syms x;

% 定义多项式
p_n5 = 1 - 1.110223e-016*x - 3.542983*x^2 + 4.440892e-016*x^3 + 2.746498*x^4;
p_n10 = 0.7308217 - 7.392743e-016*x - 4.811625*x^2 + 5.095772e-015*x^3 + 12.61929*x^4 - 2.700829e-014*x^5 - 14.00244*x^6 + 3.820127e-014*x^7 + 5.512772*x^8 - 1.573687e-014*x^9;
p_n15 = 1 + 2.553513e-015*x - 17.36407*x^2 - 2.664535e-014*x^3 + 149.0269*x^4 + 2.131628e-013*x^5 - 646.8639*x^6 + 1.080025e-012*x^7 + 1510.606*x^8 + 1.818989e-012*x^9 - 1927.183*x^10 + 1264.416*x^12 + 2.842171e-013*x^13 - 333.619*x^14;
p_n20 = 0.9624097 - 2.545115e-015*x - 16.54218*x^2 - 5.953458e-015*x^3 + 165.4582*x^4 + 7.125567e-013*x^5 - 960.8247*x^6 + 1.237559e-011*x^7 + 3379.017*x^8 - 2.322308e-011*x^9 - 7413.453*x^10 + 3.365297e-011*x^11 + 10195.47*x^12 - 2.58698e-011*x^13 - 8534.894*x^14 + 4.561769e-011*x^15 + 3973.165*x^16 - 1.688679e-011*x^17 - 788.3263*x^18 + 3.307113e-012*x^19;

% 定义函数 f(x)
f = 1 / (1 + 25*x^2);

% 创建绘图数据
x_vals = linspace(-1, 1, 1000);
p_n5_vals = double(subs(p_n5, x, x_vals));
p_n10_vals = double(subs(p_n10, x, x_vals));
p_n15_vals = double(subs(p_n15, x, x_vals));
p_n20_vals = double(subs(p_n20, x, x_vals));
f_vals = double(subs(f, x, x_vals));

% 绘图
figure;
hold on;
plot(x_vals, p_n5_vals, '--r', 'LineWidth', 1.5, 'DisplayName', 'p(x) for n=5');
plot(x_vals, p_n10_vals, '--g', 'LineWidth', 1.5, 'DisplayName', 'p(x) for n=10');
plot(x_vals, p_n15_vals, '--b', 'LineWidth', 1.5, 'DisplayName', 'p(x) for n=15');
plot(x_vals, p_n20_vals, '--m', 'LineWidth', 1.5, 'DisplayName', 'p(x) for n=20');
plot(x_vals, f_vals, 'k', 'LineWidth', 1.5, 'DisplayName', 'f(x) = 1 / (1 + 25*x^2)');

% 添加标签和图例
xlabel('x');
ylabel('p(x) and f(x)');
legend('show');
grid on;
hold off;
