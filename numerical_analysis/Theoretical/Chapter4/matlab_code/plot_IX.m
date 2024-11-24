% Define x range from a small positive value to 1
x = linspace(0.01, 1, 100);

% Define cond_f(x)
cond_f = x ./ (exp(x) - 1);

% Define the upper bound of cond_A(x)
upper_bound = exp(x) ./ x;

% Plot cond_f(x)
figure;
subplot(2, 1, 1); % Create a subplot for cond_f(x)
plot(x, cond_f, 'b-', 'LineWidth', 1.5);
xlabel('x');
ylabel('cond_f(x)');
title('Plot of cond_f(x) on [0, 1]');
grid on;

% Plot the upper bound of cond_A(x)
subplot(2, 1, 2); % Create a subplot for the upper bound
plot(x, upper_bound, 'r--', 'LineWidth', 1.5);
xlabel('x');
ylabel('Upper bound of cond_A(x)');
title('Estimated Upper Bound of cond_A(x) on [0, 1]');
grid on;
