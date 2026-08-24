% SunPathCGF.m
% Cartesian sun-path chart for Grand Forks, ND (local solar time)
clear; clc; close all;

% Location
lat = 47.93;   % latitude [deg]

% Selected daily sun paths
dayNums   = [52 113 172];                 % Feb 21, Apr 23, Jun 21
dayLabels = {'Feb 21','Apr 23','Jun 21'};

% Time and day grids
tSolar   = linspace(0,24,1000);           % local solar time [h]
nGrid    = 1:365;                         % day of year for hourly contours
hourList = 5:19;                          % local solar-time hour contours

% Figure
figure('Color','w','Position',[100 100 760 560]);
hold on; grid on; box on;

% Daily sun-path contours (blue)
for i = 1:numel(dayNums)
    n = dayNums(i);
    [az, el] = solar_position_correct(lat,n,tSolar);
    plot(az,el,'b-','LineWidth',1.8,'HandleVisibility','off');

    % Label near solar noon
    [azNoon, elNoon] = solar_position_correct(lat,n,12);
    text(azNoon,elNoon-1.8,dayLabels{i}, ...
        'Color','b','FontSize',8,'FontWeight','bold','HorizontalAlignment','center');
end

% Hourly contours (red)
for h = hourList
    [azH, elH] = solar_position_correct(lat,nGrid,h*ones(size(nGrid)));
    plot(azH,elH,'r-','LineWidth',0.8,'HandleVisibility','off');

    % Hour labels based on Jun 21 position
    [azLab, elLab] = solar_position_correct(lat,172,h);
    if ~isnan(azLab) && ~isnan(elLab)
        if h == 12
            yOffset = 1.8;
        else
            yOffset = 1.0;
        end
        text(azLab,elLab+yOffset,sprintf('%02d:00',h), ...
            'Color','r','FontSize',8,'HorizontalAlignment','center');
    end
end

% Legend handles
hDailyLegend  = plot(nan,nan,'b-','LineWidth',1.8);
hHourlyLegend = plot(nan,nan,'r-','LineWidth',0.8);
legend([hDailyLegend hHourlyLegend],{'Daily contours','Hourly contours'}, ...
    'Location','south','Box','on');

% Axes and labels
xlabel('Solar Azimuth (degrees from north, clockwise)');
ylabel('Solar Elevation (degrees)');
title({'Cartesian Sun-Path Chart for Grand Forks, ND', ...
       'Lat = 47.93^\circ N   (Local Solar Time)'});
xlim([50 310]);
ylim([0 70]);
set(gca,'FontSize',10);

% Export figure
exportgraphics(gcf,'fig_sunpath_grand_forks_corrected.png','Resolution',300);

% Local function
function [az, el] = solar_position_correct(phi,n,tSolar)
% phi    = latitude [deg]
% n      = day number (1-365)
% tSolar = local solar time [h]
% Azimuth convention: 0=north, 90=east, 180=south, 270=west

    % Solar declination [deg]
    delta = 23.45 .* sind((360/365).*(n - 81));

    % Hour angle [deg] (negative in morning)
    H = 15 .* (tSolar - 12);

    % Solar elevation
    sinEl = sind(phi).*sind(delta) + cosd(phi).*cosd(delta).*cosd(H);
    el = asind(sinEl);

    % Prevent division by zero at low elevation
    cosEl = cosd(el);

    % Solar azimuth from south-based angle gamma
    sinGamma = cosd(delta).*sind(H) ./ cosEl;
    cosGamma = (sind(el).*sind(phi) - sind(delta)) ./ (cosEl.*cosd(phi));
    gamma = atan2d(sinGamma,cosGamma);

    % Convert to azimuth from north, clockwise
    az = mod(180 + gamma,360);

    % Mask nighttime values
    az(el <= 0) = NaN;
    el(el <= 0) = NaN;
end
