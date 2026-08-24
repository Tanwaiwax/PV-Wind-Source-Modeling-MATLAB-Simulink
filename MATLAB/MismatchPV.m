% MismatchPV.m
% Two-module series-connected PV string under uniform and non-uniform irradiance
clear; clc; close all;

% Reference module parameters
Vmpp_ref = 36;       % V
Impp_ref = 5.0;      % A
Pmpp_ref = Vmpp_ref * Impp_ref;
Voc_ref  = 45;       % V
G_full = 1000;       % W/m^2

% Reduced-irradiance module target
Pmpp_reduced_target = 104;          % W
irradianceFraction = Pmpp_reduced_target / Pmpp_ref;
G_reduced = irradianceFraction * G_full;

% Construct module I-V characteristic I(V) = Isc * [1 - (V/Voc)^m]
x_mpp = Vmpp_ref / Voc_ref;
shapeFunction = @(m) (1/(m + 1))^(1/m) - x_mpp;
mShape = fzero(shapeFunction,[1 50]);

% Impp/Isc = m/(m+1)
Impp_over_Isc = mShape/(mShape + 1);
Isc_full = Impp_ref / Impp_over_Isc;
Isc_reduced = irradianceFraction * Isc_full;

% Independent module curves
V_module = linspace(0,Voc_ref,10001);
I_module_full = max(Isc_full .* (1 - (V_module./Voc_ref).^mShape),0);
P_module_full = V_module .* I_module_full;
I_module_reduced = max(Isc_reduced .* (1 - (V_module./Voc_ref).^mShape),0);
P_module_reduced = V_module .* I_module_reduced;

% Independent-module MPPs and benchmark
[Pmpp_full,idx_full] = max(P_module_full);
Vmpp_full = V_module(idx_full);
Impp_full = I_module_full(idx_full);
[Pmpp_reduced,idx_red] = max(P_module_reduced);
Vmpp_reduced = V_module(idx_red);
Impp_reduced = I_module_reduced(idx_red);
P_independent = Pmpp_full + Pmpp_reduced;

% Uniform two-module series string (both modules at full irradiance)
I_uniform = linspace(0,Isc_full,12001);
Vmod_uniform = Voc_ref .* max(1 - I_uniform./Isc_full,0).^(1/mShape);
V_string_uniform = 2 .* Vmod_uniform;
P_string_uniform = V_string_uniform .* I_uniform;
[Pmpp_uniform,idx_uniform] = max(P_string_uniform);
Vmpp_uniform = V_string_uniform(idx_uniform);
Impp_uniform = I_uniform(idx_uniform);
[V_uniform_plot,orderU] = sort(V_string_uniform);
I_uniform_plot = I_uniform(orderU);
P_uniform_plot = P_string_uniform(orderU);

% Non-uniform two-module series string (module1 full, module2 reduced)
I_series = linspace(0,Isc_reduced,12001);
V_full_series = Voc_ref .* max(1 - I_series./Isc_full,0).^(1/mShape);
V_reduced_series = Voc_ref .* max(1 - I_series./Isc_reduced,0).^(1/mShape);
V_string_nonuniform_core = V_full_series + V_reduced_series;
P_string_nonuniform_core = V_string_nonuniform_core .* I_series;

% Extend current-limited portion toward V = 0 (plateau)
V_transition = V_string_nonuniform_core(end);
V_plateau = linspace(0,V_transition,1500);
I_plateau = Isc_reduced .* ones(size(V_plateau));
P_plateau = V_plateau .* I_plateau;

% Sort core and combine with plateau for plotting
[V_core_sorted,orderN] = sort(V_string_nonuniform_core);
I_core_sorted = I_series(orderN);
P_core_sorted = P_string_nonuniform_core(orderN);
V_nonuniform_plot = [V_plateau(1:end-1), V_core_sorted];
I_nonuniform_plot = [I_plateau(1:end-1), I_core_sorted];
P_nonuniform_plot = [P_plateau(1:end-1), P_core_sorted];

% Non-uniform series-string MPP
[Pmpp_nonuniform,idx_nonuniform] = max(P_nonuniform_plot);
Vmpp_nonuniform = V_nonuniform_plot(idx_nonuniform);
Impp_nonuniform = I_nonuniform_plot(idx_nonuniform);

% Mismatch loss
P_mismatch_loss = P_independent - Pmpp_nonuniform;
Mismatch_percent = 100 * P_mismatch_loss / P_independent;

% Display numerical results
fprintf('\nINDEPENDENT MODULE MPP VALUES\n');
fprintf('Full module   Vmpp=%.3f V, Impp=%.3f A, Pmpp=%.3f W\n',Vmpp_full,Impp_full,Pmpp_full);
fprintf('Reduced module Vmpp=%.3f V, Impp=%.3f A, Pmpp=%.3f W\n',Vmpp_reduced,Impp_reduced,Pmpp_reduced);
fprintf('Independent benchmark P = %.3f W\n\n',P_independent);

fprintf('UNIFORM SERIES STRING: Vmpp=%.3f V, Impp=%.3f A, Pmpp=%.3f W\n',Vmpp_uniform,Impp_uniform,Pmpp_uniform);
fprintf('NON-UNIFORM SERIES STRING: Vmpp=%.3f V, Impp=%.3f A, Pmpp=%.3f W\n\n',Vmpp_nonuniform,Impp_nonuniform,Pmpp_nonuniform);

fprintf('MISMATCH RESULTS: Independent=%.3f W, Series MPP=%.3f W, Loss=%.3f W (%.2f%%)\n', ...
    P_independent, Pmpp_nonuniform, P_mismatch_loss, Mismatch_percent);

% Rounded values for manuscript
fprintf('\nMANUSCRIPT-ROUNDED VALUES\n');
fprintf('Full module MPP        = %.0f W\n',Pmpp_full);
fprintf('Reduced module MPP     = %.0f W\n',Pmpp_reduced);
fprintf('Independent benchmark  = %.0f W\n',P_independent);
fprintf('Non-uniform string MPP = %.0f W\n',Pmpp_nonuniform);
fprintf('Mismatch loss          = %.0f W (≈ %.0f%%)\n',P_mismatch_loss,Mismatch_percent);

% Plot
figure('Color','w','Position',[100 100 760 650]);
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
plot(V_uniform_plot,I_uniform_plot,'b-','LineWidth',1.8,'DisplayName','Uniform');
hold on;
plot(V_nonuniform_plot,I_nonuniform_plot,'r--','LineWidth',1.8,'DisplayName','Non-uniform');
plot(Vmpp_uniform,Impp_uniform,'bo','MarkerFaceColor','b','MarkerSize',6,'HandleVisibility','off');
plot(Vmpp_nonuniform,Impp_nonuniform,'rs','MarkerFaceColor','r','MarkerSize',6,'HandleVisibility','off');
xlabel('String Voltage, V [V]'); ylabel('String Current, I [A]');
title('(a) Current-Voltage (I-V) Characteristics');
xlim([0 95]); ylim([0 6]); grid on; box on; legend('Location','southwest');

ax2 = nexttile;
plot(V_uniform_plot,P_uniform_plot,'b-','LineWidth',1.8,'DisplayName','Uniform');
hold on;
plot(V_nonuniform_plot,P_nonuniform_plot,'r--','LineWidth',1.8,'DisplayName','Non-uniform');
plot(Vmpp_uniform,Pmpp_uniform,'bo','MarkerFaceColor','b','MarkerSize',6,'HandleVisibility','off');
plot(Vmpp_nonuniform,Pmpp_nonuniform,'rs','MarkerFaceColor','r','MarkerSize',6,'HandleVisibility','off');
text(Vmpp_uniform-26,Pmpp_uniform+15,sprintf('Uniform MPP: %.0f W',Pmpp_uniform),'FontSize',9);
text(Vmpp_nonuniform-27,Pmpp_nonuniform-28,sprintf('Non-uniform MPP: %.0f W',Pmpp_nonuniform),'FontSize',9);
xlabel('String Voltage, V [V]'); ylabel('String Power, P [W]');
title('(b) Power-Voltage (P-V) Characteristics');
xlim([0 95]); ylim([0 400]); grid on; box on; legend('Location','northwest');

linkaxes([ax1 ax2],'x');
title(t,'Two-Module Series PV String under Uniform and Non-Uniform Irradiance','FontWeight','bold');

% Export figure
exportgraphics(gcf,'Figure9_Series_String_Mismatch.png','Resolution',600);
