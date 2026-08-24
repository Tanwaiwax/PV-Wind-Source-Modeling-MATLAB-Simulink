%% Rayleigh wind-speed PDF, cubic weighting, and power-contribution weighting

clear; clc; close all;

%% Parameters

vbar = 7.55;           % Average wind speed [m/s]
vmax = 25;             % maximum wind speed for plotting [m/s]
v = linspace(0, vmax, 1000);

%% Rayleigh scale parameter

c = 2*vbar/sqrt(pi);   % c = 8.52 m/s

fprintf('Mean wind speed, vbar = %.2f m/s\n', vbar);
fprintf('Rayleigh scale parameter, c = %.2f m/s\n', c);

%% Rayleigh probability density function

f_rayleigh = (2*v./c.^2).*exp(-(v./c).^2);

%% Normalized cubic wind-power weighting

cubic_norm = v.^3 ./ max(v.^3);

%% Normalized Rayleigh-weighted wind-power contribution

power_contribution = (v.^3).*f_rayleigh;
power_contribution_norm = ...
    power_contribution ./ max(power_contribution);

%% Numerical check of Rayleigh third-moment relation

v3_avg_num = trapz(v, (v.^3).*f_rayleigh);
ratio_num = v3_avg_num/(vbar^3);

fprintf('Numerical E[v^3]/vbar^3 = %.3f\n', ratio_num);
fprintf('Expected Rayleigh value = 1.91\n');

%% Plot

figure('Color','w','Position',[100 100 950 500]);

yyaxis left
plot(v, f_rayleigh, 'b-', 'LineWidth', 2.0);
ylabel('Rayleigh probability density, f(v)');
ylim([0 1.15*max(f_rayleigh)]);

yyaxis right
plot(v, cubic_norm, 'r-', 'LineWidth', 2.0);
hold on;
plot(v, power_contribution_norm, 'k--', 'LineWidth', 2.0);
ylabel('Normalized weighting');
ylim([0 1.05]);

xlabel('Wind speed, v (m/s)');
title('Rayleigh Wind-Speed PDF, Cubic Weighting, and Power-Contribution Weighting');

grid on;
box on;
xlim([0 vmax]);

legend('Rayleigh PDF, f(v)', ...
       'Normalized cubic weighting, v^3/max(v^3)', ...
       'Normalized power contribution, v^3f(v)/max[v^3f(v)]', ...
       'Location','northeast', ...
       'Orientation','vertical');

%% Export figure

exportgraphics(gcf, ...
    'fig_rayleigh_wind_power_threecurves.png', ...
    'Resolution', 300);