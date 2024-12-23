% B_plot.m

clear;

% Define the symbolic variable
syms x;

% Define f(x) as a symbolic expression
f = 1 / (1 + x^2);

% Define interpolation polynomials
p2 = 1 - 0.03846154 * x^2;
p4 = 1 + 2.775558e-017 * x - 0.1710875 * x^2 + 3.469447e-018 * x^3 + 0.00530504 * x^4;
p6 = 1 - 5.551115e-017 * x - 0.3513637 * x^2 + 0.0335319 * x^4 - 8.673617e-019 * x^5 - 0.0008406327 * x^6;
p8 = 1 - 1.387779e-017 * x - 0.5281214 * x^2 - 9.714451e-017 * x^3 + 0.09818753 * x^4 + 8.673617e-018 * x^5 - 0.006580161 * x^6 + 0.0001374446 * x^8;

% Create a range of x values for plotting
x_vals = linspace(-5, 5, 1000);

% Convert symbolic expressions to function handles for numerical evaluation
f_func = matlabFunction(f);
p2_func = matlabFunction(p2);
p4_func = matlabFunction(p4);
p6_func = matlabFunction(p6);
p8_func = matlabFunction(p8);

% Calculate values for plotting
f_vals = f_func(x_vals);
p2_vals = p2_func(x_vals);
p4_vals = p4_func(x_vals);
p6_vals = p6_func(x_vals);
p8_vals = p8_func(x_vals);

% Plot the results
figure;
plot(x_vals, f_vals, 'k', 'LineWidth', 2); % Original function f(x)
hold on;
plot(x_vals, p2_vals, '--r', 'LineWidth', 1.5); % n=2 polynomial
plot(x_vals, p4_vals, '--g', 'LineWidth', 1.5); % n=4 polynomial
plot(x_vals, p6_vals, '--b', 'LineWidth', 1.5); % n=6 polynomial
plot(x_vals, p8_vals, '--m', 'LineWidth', 1.5); % n=8 polynomial

% Add legend and labels
legend('f(x) = 1 / (1 + x^2)', 'p(x), n=2', 'p(x), n=4', 'p(x), n=6', 'p(x), n=8', 'Location', 'Best');
xlabel('x');
ylabel('y');
grid on;

