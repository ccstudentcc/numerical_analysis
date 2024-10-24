% F_plot.m

% Clear workspace and close figures
clear; close all;

% Read and plot data for m = 10
points_10 = dlmread('..\data\data_10.txt', ',', 0, 0);
figure;
scatter(points_10(:, 1), points_10(:, 2), 'filled');
xlabel('X-axis');
ylabel('Y-axis');

grid on;

% Read and plot data for m = 40
points_40 = dlmread('..\data\data_40.txt', ',', 0, 0);
figure;
scatter(points_40(:, 1), points_40(:, 2), 'filled');
xlabel('X-axis');
ylabel('Y-axis');

grid on;

% Read and plot data for m = 160
points_160 = dlmread('..\data\data_160.txt', ',', 0, 0);
figure;
scatter(points_160(:, 1), points_160(:, 2), 'filled');
xlabel('X-axis');
ylabel('Y-axis');

grid on;


