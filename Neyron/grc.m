close all;
clc;
clear;

FontSize = 10;
filename1 = '/users/q/source/repos/Neyron/Neyron/results.txt';
Data = importdata( [ filename1 ] );
time1   = Data(:, 1);
V1   = Data(:, 2);
time2    = Data(:, 3);
V2      = Data(:, 4);

fig_Ca_IP3_z = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
subplot(1, 1, 1);
plot(time1, V1, 'r', 'LineWidth', 1);

hold on;
grid on;
ylim([0 0,5]);
xlim([0 180.0]);
set(gca, 'YTick', [0: 0.1: 0.7]);
