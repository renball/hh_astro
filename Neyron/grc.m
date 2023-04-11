close all;
clc;
clear;

FontSize = 10;
filename1 = '/users/q/source/repos/Neyron/Neyron/results.txt';
Data = importdata( [ filename1 ] );
time   = Data(:, 1);
Ca     = Data(:, 2);
IP3    = Data(:, 3);
z      = Data(:, 4);

fig_Ca_IP3_z = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
subplot(1, 1, 1);
plot(time, Ca, 'r', 'LineWidth', 1);
hold on;
grid on;
ylim([0 0.7]);
xlim([0 1.0]);
set(gca, 'YTick', [0: 0.1: 0.7]);
