% SingleDiode_IV_PV_Analytical.m
% Analytical I-V and P-V characteristics of a representative single-diode PV cell
clear; clc; close all;

% Parameters
Isc = 1.2;          % short-circuit current [A]
I0  = 2.0e-9;       % reverse saturation current [A]
A   = 40;           % diode exponential constant [1/V]
Rp  = 12;           % shunt resistance [ohm]

% Voltage sweep
V = linspace(0,0.60,1000);   % voltage [V]

% Ideal single-diode model (no shunt)
I_ideal = Isc - I0 .* (exp(A .* V) - 1);

% Single-diode with finite shunt resistance Rp
I_Rp = Isc - I0 .* (exp(A .* V) - 1) - V ./ Rp;

% Remove non-physical negative current after open-circuit
I_ideal(I_ideal < 0) = 0;
I_Rp(I_Rp < 0)       = 0;

% Power curves
P_ideal = V .* I_ideal;
P_Rp    = V .* I_Rp;

% Maximum power points
[Pmax_ideal, idx_ideal] = max(P_ideal);
Vmpp_ideal = V(idx_ideal);
Impp_ideal = I_ideal(idx_ideal);

[Pmax_Rp, idx_Rp] = max(P_Rp);
Vmpp_Rp = V(idx_Rp);
Impp_Rp = I_Rp(idx_Rp);

% Display results
fprintf('Ideal model MPP:\n');
fprintf('Vmpp = %.4f V, Impp = %.4f A, Pmpp = %.4f W\n\n', Vmpp_ideal, Impp_ideal, Pmax_ideal);

fprintf('Finite Rp model MPP:\n');
fprintf('Vmpp = %.4f V, Impp = %.4f A, Pmpp = %.4f W\n\n', Vmpp_Rp, Impp_Rp, Pmax_Rp);

% Plot I-V and P-V with two y-axes
figure('Color','w','Position',[100 100 900 430]);

yyaxis left
plot(V, I_ideal, '-', 'LineWidth', 1.6); hold on;
plot(V, I_Rp,    '--', 'LineWidth', 1.6);
ylabel('Current I [A]');
ylim([0 1.25]);
ax = gca;
ax.YColor = [0 0.4470 0.7410];

yyaxis right
plot(V, P_ideal, '-', 'LineWidth', 1.6);
plot(V, P_Rp,    '--', 'LineWidth', 1.6);
ylabel('Power P [W]');
ylim([0 0.52]);
ax.YColor = [0.8500 0.3250 0.0980];

xlabel('Voltage V [V]');
title('Analytical I-V and P-V Characteristics of Single-Diode PV Cell');

grid on; box on; xlim([0 0.60]);

legend('I-V ideal', 'I-V with R_p', 'P-V ideal', 'P-V with R_p', 'Location', 'south');

% Export figure
exportgraphics(gcf, 'fig_single_diode_iv_pv.png', 'Resolution', 300);
