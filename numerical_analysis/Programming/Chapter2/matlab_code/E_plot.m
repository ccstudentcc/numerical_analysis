% B_plot_symbolic.m

clear;

% Define the symbolic variable
syms x;


% Define interpolation polynomials
sp1 = 4.1477e-05*x^6 + -0.00371557*x^5 + 0.128281*x^4 + -2.11512*x^3 + 16.2855*x^2 + -43.0127*x^1 + 6.67*x^0;
sp2 = 6.67-5.85018*x+2.982271*x^2-0.4242825*x^3+0.02658578*x^4-0.0007774732*x^5+8.676802e-006*x^6;

% Create a range of x values for plotting
x_vals = linspace(0, 30, 1000);

% Convert symbolic expressions to function handles for numerical evaluation
sp1_func = matlabFunction(sp1);
sp2_func = matlabFunction(sp2);


% Calculate values for plotting
sp1_vals = sp1_func(x_vals);
sp2_vals = sp2_func(x_vals);


% Plot the results
figure;
hold on;
plot(x_vals, sp1_vals, 'r', 'LineWidth', 1.5); 
plot(x_vals, sp2_vals, 'g', 'LineWidth', 1.5); 


% Add legend and labels
legend( 'sp1', 'sp2', 'Location', 'Best');
xlabel('day');
ylabel('weight');
grid on;

