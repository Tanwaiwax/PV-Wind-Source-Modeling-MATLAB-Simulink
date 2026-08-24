% NSRDBGFPV.m
% High-irradiance comparison: NSRDB-derived vs analytical clear-sky POA
clear; clc; close all;

% Input file
filename = 'NSRDB PVdata HESS.csv';

% Site and collector
lat  = 47.93; lon = -97.02; tz = -6; elev = 254;
tilt = lat; surfAz = 180; rho_g = 0.20;

% Selection thresholds
minSolarAltitude = 5;
poaMinThreshold = 850;

% Read data (first two rows metadata)
T = readtable(filename,'NumHeaderLines',2,'VariableNamingRule','preserve');
yearCol   = T.Year; monthCol  = T.Month; dayCol    = T.Day;
hourCol   = T.Hour; minuteCol = T.Minute;
DHI = T.DHI; DNI = T.DNI; GHI = T.GHI;

% Day-of-year using a fixed geometry year
geometryYear = 2018;
dt = datetime(geometryYear, monthCol, dayCol, hourCol, minuteCol, 0);
dayNumber = day(dt,'dayofyear');
localTime = hourCol + minuteCol/60;

% Solar-time correction
B = 360*(dayNumber - 81)/364;
EoT = 9.87*sind(2*B) - 7.53*cosd(B) - 1.5*sind(B);
LSTM = 15*tz;
TC = 4*(lon - LSTM) + EoT;
solarTime = localTime + TC/60;

% Solar geometry
delta = 23.45*sind((360/365).*(dayNumber - 81));
H = 15*(solarTime - 12);
sinBeta = sind(lat).*sind(delta) + cosd(lat).*cosd(delta).*cosd(H);
sinBeta = max(min(sinBeta,1),-1);
beta = asind(sinBeta);
cosBeta = cosd(beta);
cosBeta(abs(cosBeta) < 1e-8) = NaN;
sinGamma = cosd(delta).*sind(H)./cosBeta;
cosGamma = (sind(beta).*sind(lat) - sind(delta)) ./ (cosBeta.*cosd(lat));
gammaSouth = atan2d(sinGamma,cosGamma);
solarAz = mod(180 + gammaSouth,360);

% Incidence angle on tilted collector
cosTheta = sind(beta).*cosd(tilt) + cosd(beta).*sind(tilt).*cosd(solarAz - surfAz);
cosTheta = max(cosTheta,0);

% NSRDB-derived POA
Ibc_data = DNI .* cosTheta;
Idc_data = DHI .* ((1 + cosd(tilt))/2);
Irc_data = GHI .* rho_g .* ((1 - cosd(tilt))/2);
Gpoa_data = Ibc_data + Idc_data + Irc_data;
Gpoa_data(beta <= 0) = 0;

% Analytical clear-sky POA (Masters-style)
As = 1160 + 75*sind((360/365).*(dayNumber - 275));
k  = 0.174 + 0.035*sind((360/365).*(dayNumber - 100));
C  = 0.095 + 0.04*sind((360/365).*(dayNumber - 100));
IB_clear = zeros(size(beta));
dayIdx = beta > 0;
IB_clear(dayIdx) = As(dayIdx).*exp(-k(dayIdx)./sind(beta(dayIdx)));
IG_clear = IB_clear .* (sind(beta) + C);
Ibc_clear = IB_clear .* cosTheta;
Idc_clear = C .* IB_clear .* ((1 + cosd(tilt))/2);
Irc_clear = IG_clear .* rho_g .* ((1 - cosd(tilt))/2);
Gpoa_clear = Ibc_clear + Idc_clear + Irc_clear;
Gpoa_clear(beta <= 0) = 0;

% Select high-irradiance daytime samples (selection uses only NSRDB-derived POA)
idx = beta > minSolarAltitude & Gpoa_data >= poaMinThreshold & isfinite(Gpoa_data) & isfinite(Gpoa_clear);
x = Gpoa_data(idx); % NSRDB-derived POA
y = Gpoa_clear(idx); % Analytical clear-sky POA
if isempty(x), error('No samples satisfy the selected irradiance criteria.'); end

% Error statistics: analytical - NSRDB
err = y - x;
npts = numel(err);
MBE = mean(err);
RMSE = sqrt(mean(err.^2));
MAE = mean(abs(err));
R = corr(x,y);

% Axis limits
allVals = [x(:); y(:)];
axisMin = floor(min(allVals)/50)*50;
axisMax = ceil(max(allVals)/50)*50;
axisMin = max(0,axisMin - 50);
axisMax = axisMax + 50;

% Plot
figure('Color','w','Position',[100 100 760 620]);
scatter(x,y,18,'MarkerFaceColor',[0.20 0.20 0.20],'MarkerEdgeColor',[0.20 0.20 0.20], ...
    'MarkerFaceAlpha',0.35,'MarkerEdgeAlpha',0.20); hold on;
plot([axisMin axisMax],[axisMin axisMax],'r--','LineWidth',1.8);
xlabel('NSRDB-derived POA irradiance [W/m^2]','FontSize',11);
ylabel('Analytical clear-sky POA irradiance [W/m^2]','FontSize',11);
title({'Plane-of-Array (POA) Irradiance Comparison','High-Irradiance Daytime Samples, Grand Forks, ND'}, ...
    'FontSize',12,'FontWeight','bold');
xlim([axisMin axisMax]); ylim([axisMin axisMax]);
axis square; grid on; box on; set(gca,'FontSize',10);
legend('POA sample pairs','1:1 reference line','Location','northwest');

% Statistics annotation
statsText = sprintf(['NSRDB POA \\geq %.0f W/m^2\nn = %d\nMBE = %.1f W/m^2\nRMSE = %.1f W/m^2'], ...
    poaMinThreshold, npts, MBE, RMSE);
text(0.97,0.05,statsText,'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','bottom','FontSize',9,'BackgroundColor','w','EdgeColor',[0.45 0.45 0.45],'Margin',6);

% Export
exportgraphics(gcf,'fig_poa_high_irradiance_comparison.png','Resolution',600);
