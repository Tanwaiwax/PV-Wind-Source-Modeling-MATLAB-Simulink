clear; clc; close all;

%% Synthetic POA datasets
rng(10);

N = 6000;

% Dataset 1 centered near 950 W/m^2
x = 950 + 8*randn(N,1);

% Dataset 2 correlated with dataset 1
y = x + 4*randn(N,1);

% Constrain both to 900–1000 W/m^2
x = max(min(x,1000),900);
y = max(min(y,1000),900);

%% Plot
figure('Color','w','Position',[100 100 650 520]);

% Black particle cloud
scatter(x,y,6,'k','filled', ...
    'MarkerFaceAlpha',0.08, ...
    'MarkerEdgeAlpha',0.08);

hold on;

% 1:1 reference line
plot([900 1000],[900 1000],'r--','LineWidth',1.8);

xlabel('Synthetic Plane-of-Array (POA) irradiance set 1 [W/m^2]', ...
       'FontSize',12);

ylabel('Synthetic Plane-of-Array (POA) irradiance set 2 [W/m^2]', ...
       'FontSize',12);

title('Synthetic Plane-of-Array (POA) Irradiance Comparison', ...
      'FontSize',13,'FontWeight','bold');

xlim([900 1000]);
ylim([900 1000]);

grid on;
box on;
set(gca,'FontSize',11);

%% Export
exportgraphics(gcf,'fig_poa_comparison_corrected.png','Resolution',300);