% Define the knots t_i
%t = 0:5;  % Example knots with a range of values
i=0;
x = linspace(i-3, i+4, 1000);  % A fine grid for plotting

% Define the quadratic B-spline basis function B^2_i(x)
B_i2 = zeros(size(x));
for j = 1:length(x)
    if x(j) > (i-1) && x(j) <= (i)
        B_i2(j) = ((x(j) - (i-1))^2) / (((i) - (i-1)) * ((i+1) - (i-1)));
    elseif x(j) > (i) && x(j) <= (i+1)
        B_i2(j) = (x(j) - (i-1))*((i+1)-x(j)) / (((i+1) - (i-1)) * ((i+1) - (i))) + (x(j) - (i))*((i+2)-x(j)) / (((i+2) - (i)) * ((i+1) - (i)));
    elseif x(j) > (i+1) && x(j) <= (i+2)
        B_i2(j) = (((i+2) - x(j))^2) / (((i+2) - (i+1)) * ((i+2) - (i)));
    else
        B_i2(j) = 0; % Outside the support
    end
end

% for j = 1:length(x)
%     if x(j) > t(i-1) && x(j) <= t(i)
%         B_i2(j) = ((x(j) - t(i-1))^2) / ((t(i) - t(i-1)) * (t(i+1) - t(i-1)));
%     elseif x(j) > t(i) && x(j) <= t(i+1)
%         B_i2(j) = (x(j) - t(i-1))*(t(i+1)-x(j)) / ((t(i+1) - t(i-1)) * (t(i+1) - t(i))) - (x(j) - t(i))*(t(i+2)-x(j)) / ((t(i+2) - t(i)) * (t(i+1) - t(i)));
%     elseif x(j) > t(i+1) && x(j) <= t(i+2)
%         B_i2(j) = ((t(i+2) - x(j))^2) / ((t(i+2) - t(i+1)) * (t(i+2) - t(i)));
%     else
%         B_i2(j) = 0; % Outside the support
%     end
% end

% Plot the B-spline basis function
plot(x, B_i2, 'LineWidth', 2);
xlabel('x');
ylabel('B^2_i(x)');
grid on;

% Adjust x-axis to display i-3, i-2, i-1, i, i+1, i+2, i+3
xlim([i-3, i+4]);  % Set x-axis limits to show the relevant range
xticks(i-3:i+4);  % Set ticks at i-3, i-2, ..., i+3
xticklabels(arrayfun(@(n) sprintf('i%+d', n - i), i+1:i+3, 'UniformOutput', false)); % Custom labels i-3, i-2, ..., i+3
xticklabels({'i-3','i-2','i-1','i','i+1','i+2','i+3','i+4'});