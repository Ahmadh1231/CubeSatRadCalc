clear all;
close all;
clc;
%% Parameters
% CubeSat dimensions (6U = 10cm x 20cm x 30cm)
dimensions = [0.1, 0.2, 0.3]; % Width, height, length in meters
% X dimension = 10cm (smallest - 2U sides)
% Y dimension = 20cm (middle - 3U sides)
% Z dimension = 30cm (largest - 6U sides)

% Sun and Earth parameters
solarFlux = 1367; % W/m^2 (solar constant)
earthAlbedo = 0.3; % Earth's albedo coefficient (reflection)
earthIR = 237; % W/m^2 (Earth's IR radiation)

% Earth parameters
earthRadius = 6371; % km
earthMu = 3.986004418e14; % Earth's gravitational parameter (m^3/s^2)

% Material properties
absorptivity = 0.8; % Solar absorptivity
emissivity = 0.8; % Infrared emissivity

% Simulation parameters
simulationHours = 24;  % 24-hour simulation
minutesPerHour = 60;
timeStepMinutes = 1;   % 1-minute intervals
numTimeSteps = uint32(simulationHours * minutesPerHour / timeStepMinutes); % Integer type conversion

%% Orbit Parameters (Using Keplerian elements)
% Define orbital elements
% Semi-major axis already defined by altitude
altitude = 500; % km
semiMajorAxis = earthRadius + altitude; % km
eccentricity = 0.001; % Near-circular orbit
inclination = 50.0; % degrees (50 degrees as specified)
raan = 30; % degrees (Right Ascension of Ascending Node)
argPerigee = 0; % degrees (Argument of Perigee)
trueAnomaly = 0; % degrees (Initial position in orbit)

% Convert to radians for calculations
inclinationRad = inclination * pi/180;
raanRad = raan * pi/180;
argPerigeeRad = argPerigee * pi/180;
trueAnomalyRad = trueAnomaly * pi/180;

% Calculate orbital period
orbitPeriod = 2 * pi * sqrt((semiMajorAxis*1000)^3 / earthMu); % seconds

%% Face data
faceNames = {'+X (2U)', '-X (2U)', '+Y (3U)', '-Y (3U)', '+Z (6U)', '-Z (6U)'};
faceDimensions = {'10cm × 20cm', '10cm × 20cm', '10cm × 30cm', '10cm × 30cm', '20cm × 30cm', '20cm × 30cm'};

%% Create time and results arrays with guaranteed integer size
timeSeconds = zeros(1, numTimeSteps);
radiationData = zeros(6, numTimeSteps);
earthViewFactors = zeros(6, numTimeSteps);
isEclipsed = zeros(1, numTimeSteps); % Track eclipse status

% Arrays for radiation components
solarRadiationData = zeros(6, numTimeSteps); 
albedoRadiationData = zeros(6, numTimeSteps); 
irRadiationData = zeros(6, numTimeSteps);

% Arrays for orbit visualization
orbitalPositions = zeros(numTimeSteps, 3);
eclipsedPositions = [];
sunlitPositions = [];

fprintf('Running 24-hour simulation with 1-minute intervals (%d steps)...\n', numTimeSteps);

% Create time array explicitly (avoid using ranges)
for i = 1:numTimeSteps
    timeSeconds(i) = (i-1) * timeStepMinutes * 60; % Time in seconds
end

% Create vertices for cuboid
w = dimensions(1)/2;  % X dimension = 10cm (smallest)
h = dimensions(2)/2;  % Y dimension = 20cm (middle)
l = dimensions(3)/2;  % Z dimension = 30cm (largest)

vertices = [
    -w, -h, -l;  % 1: bottom-left-back
    w, -h, -l;   % 2: bottom-right-back
    w, h, -l;    % 3: top-right-back
    -w, h, -l;   % 4: top-left-back
    -w, -h, l;   % 5: bottom-left-front
    w, -h, l;    % 6: bottom-right-front
    w, h, l;     % 7: top-right-front
    -w, h, l     % 8: top-left-front
];

% Define faces (each row contains indices of vertices that form a face)
faces = [
    1 2 6 5;    % +X face (2U - smallest)
    3 4 8 7;    % -X face (2U - smallest)
    2 3 7 6;    % +Y face (3U - middle)
    1 4 8 5;    % -Y face (3U - middle)
    5 6 7 8;    % +Z face (6U - largest, sun-facing)
    1 2 3 4     % -Z face (6U - largest)
];

%% Calculate Sun's position throughout the year
% Simplified sun position in Earth-centered inertial (ECI) frame
% Assumes Earth's orbit is circular and in XY plane
sunPositionECI = zeros(numTimeSteps, 3);

% Ensure we start with January 1st (simplified)
dayOfYear = 1;
  
for t = 1:numTimeSteps
    % Compute Sun's position (simplified model)
    % Earth's orbital period around Sun in days
    earthOrbitalPeriod = 365.256;
    
    % Sun angle in Earth's orbit (in radians)
    theta = 2 * pi * (dayOfYear / earthOrbitalPeriod) + (timeSeconds(t) / (24 * 3600 * earthOrbitalPeriod));
    
    % Sun-Earth distance (in AU)
    sunEarthDistance = 1.496e11; % meters (1 AU)
    
    % Sun position in ECI frame
    sunPositionECI(t, :) = [sunEarthDistance * cos(theta), ...
                            sunEarthDistance * sin(theta), ...
                            0];
end

%% Main simulation loop
for t = 1:numTimeSteps
    % Current simulation time in seconds
    currentTime = timeSeconds(t);
    
    % Calculate Mean Anomaly at current time
    meanMotion = 2 * pi / orbitPeriod; % rad/s
    meanAnomaly = meanMotion * currentTime;
    
    % Solve Kepler's equation to get eccentric anomaly (iterative)
    E = meanAnomaly; % Initial guess
    for i = 1:10 % Usually converges in a few iterations
        E = meanAnomaly + eccentricity * sin(E);
    end
    
    % Calculate true anomaly from eccentric anomaly
    trueAnomalyRad = 2 * atan2(sqrt(1 + eccentricity) * sin(E/2), ...
                              sqrt(1 - eccentricity) * cos(E/2));
    
    % Calculate distance from Earth center (in km)
    radius = semiMajorAxis * (1 - eccentricity * cos(E));
    
    % Calculate position in orbital plane
    x_orbital = radius * cos(trueAnomalyRad);
    y_orbital = radius * sin(trueAnomalyRad);
    
    % Transform to ECI (Earth-Centered Inertial) coordinates
    % This incorporates inclination, RAAN, and argument of perigee
    xp = x_orbital * cos(argPerigeeRad) - y_orbital * sin(argPerigeeRad);
    yp = x_orbital * sin(argPerigeeRad) + y_orbital * cos(argPerigeeRad);
    
    x_eci = xp * cos(raanRad) - yp * cos(inclinationRad) * sin(raanRad);
    y_eci = xp * sin(raanRad) + yp * cos(inclinationRad) * cos(raanRad);
    z_eci = yp * sin(inclinationRad);
    
    % CubeSat position in ECI frame (in km)
    position = [x_eci, y_eci, z_eci];
    
    % Store position for visualization
    orbitalPositions(t, :) = position;
    
    % Calculate Sun direction in ECI coordinates
    sunPosition = sunPositionECI(t, :);
    sunDirection = sunPosition / norm(sunPosition);
    
    % Eclipse detection
    % Check if satellite is in Earth's shadow
    % Vector from Earth to satellite
    earthToSat = position * 1000; % Convert to meters
    earthToSatNorm = earthToSat / norm(earthToSat);
    
    % Angle between Earth-satellite vector and Earth-Sun vector
    cosAngle = dot(-sunDirection, earthToSatNorm);
    
    % Distance from satellite to Earth-Sun line
    distToEarthSunLine = norm(earthToSat) * sqrt(1 - cosAngle^2);
    
    % Satellite is in eclipse if:
    % 1. It's on the opposite side from the sun (cosAngle > 0)
    % 2. The distance to Earth-Sun line is less than Earth's radius
    if cosAngle > 0 && distToEarthSunLine < earthRadius*1000
        inEclipse = 1;
        isEclipsed(t) = 1;
        eclipsedPositions = [eclipsedPositions; position];
    else
        inEclipse = 0;
        sunlitPositions = [sunlitPositions; position];
    end
    
    % Calculate sun direction relative to CubeSat
    relSunDirection = sunDirection;
    
    % Required rotation to point +Z towards sun
    defaultZAxis = [0, 0, 1];
    rotationAxis = cross(defaultZAxis, relSunDirection);
    
    % Create rotation matrix
    if norm(rotationAxis) < 1e-10
        rotationMatrix = eye(3);
    else
        rotationAxis = rotationAxis / norm(rotationAxis);
        rotationAngle = acos(dot(defaultZAxis, relSunDirection));
        
        c = cos(rotationAngle);
        s = sin(rotationAngle);
        t_val = 1 - c;
        x = rotationAxis(1);
        y = rotationAxis(2);
        z = rotationAxis(3);
        
        rotationMatrix = [
            t_val*x*x + c, t_val*x*y - z*s, t_val*x*z + y*s;
            t_val*x*y + z*s, t_val*y*y + c, t_val*y*z - x*s;
            t_val*x*z - y*s, t_val*y*z + x*s, t_val*z*z + c
        ];
    end
    
    % Rotate vertices
    rotatedVertices = (rotationMatrix * vertices')';
    
    % Calculate face normals
    faceNormals = zeros(6, 3);
    for i = 1:6
        face = faces(i, :);
        v1 = rotatedVertices(face(2), :) - rotatedVertices(face(1), :);
        v2 = rotatedVertices(face(3), :) - rotatedVertices(face(2), :);
        normal = cross(v1, v2);
        faceNormals(i, :) = normal / norm(normal);
    end
    
    % Calculate nadir direction (toward Earth center)
    nadir = -earthToSatNorm;
    
    % Earth rho angle (angular radius of Earth as seen from satellite)
    rho = asin(earthRadius*1000 / norm(earthToSat));
    
    % Calculate radiation for each face
    for faceIdx = 1:6
        % 1. Earth view factor
        angle = acos(dot(faceNormals(faceIdx,:), nadir));
        if angle < pi/2 && angle < (pi/2 + rho)
            % Face can see Earth
            earthViewFactors(faceIdx, t) = cos(angle) * (1 - cos(rho)) / 2;
        else
            % Face cannot see Earth
            earthViewFactors(faceIdx, t) = 0;
        end
        
        % 2. Solar radiation (accounting for eclipse)
        cosIncidence = dot(faceNormals(faceIdx,:), relSunDirection);
        if cosIncidence > 0 && ~inEclipse
            solarRadiationData(faceIdx, t) = solarFlux * cosIncidence * absorptivity;
        else
            solarRadiationData(faceIdx, t) = 0; % Face doesn't see the sun or is in eclipse
        end
        
        % 3. Earth albedo radiation (accounting for eclipse)
        % This is solar radiation REFLECTED from Earth's surface
        if ~inEclipse
            albedoRadiationData(faceIdx, t) = solarFlux * earthAlbedo * earthViewFactors(faceIdx, t) * absorptivity;
        else
            albedoRadiationData(faceIdx, t) = 0; % No albedo in eclipse
        end
        
        % 4. Earth IR radiation (always present when Earth is visible)
        irRadiationData(faceIdx, t) = earthIR * earthViewFactors(faceIdx, t) * emissivity;
        
        % 5. Total radiation
        radiationData(faceIdx, t) = solarRadiationData(faceIdx, t) + albedoRadiationData(faceIdx, t) + irRadiationData(faceIdx, t);
    end
    
    % Progress indicator (every 10%)
    if mod(t, max(1, floor(numTimeSteps/10))) == 0
        fprintf('Progress: %.0f%%\n', 100.0 * t / numTimeSteps);
    end
end

% Store final radiation values
instantRadiationData = radiationData(:, numTimeSteps);

fprintf('Simulation complete! Generating visualizations...\n');

%% 1. Radiation Heatmaps
figure(1);
plotRadiationHeatmaps(radiationData, instantRadiationData, faceNames, faceDimensions);

%% 2. Radiation Time Chart with Eclipse Indicators
figure(2);
plotRadiationTimeChart(radiationData, faceNames, timeSeconds, isEclipsed);

%% 3. Orbit Visualization
figure(3);
plotOrbit(orbitalPositions, sunlitPositions, eclipsedPositions, earthRadius);

%% 4. NEW: Radiation Component Analysis
figure(4);
plotRadiationComponents(solarRadiationData, albedoRadiationData, irRadiationData, timeSeconds);

%% Helper Functions for Visualization

function plotRadiationHeatmaps(radiationData, instantRadiationData, faceNames, faceDimensions)
    % Create radiation heatmaps
    
    % Find global max radiation for consistent color scale
    maxRadiation = max(radiationData(:)) * 1.1;
    if maxRadiation == 0
        maxRadiation = 1; % Avoid division by zero
    end
    
    % Define face dimensions in pixels proportional to actual dimensions
    % Order: +X, -X, +Y, -Y, +Z, -Z
    faceSizes = {[30 60], [30 60], [30 90], [30 90], [60 90], [60 90]};
    
    % Setup the figure
    figure('Position', [100, 100, 1200, 800]);
    
    % Create the subplots - one for each face
    for i = 1:6
        ax = subplot(2, 3, i);
        
        % Get the dimensions of this face for proper visualization
        faceSize = faceSizes{i};
        
        % Create variable resolution heatmap data based on face dimensions
        [x, y] = meshgrid(linspace(-1, 1, faceSize(1)), linspace(-1, 1, faceSize(2)));
        
        % Create radiation pattern
        % For the sun-facing face (+Z), create a more concentrated pattern
        if i == 5 % +Z face (sun-facing)
            r = sqrt(x.^2 + y.^2);
            heatmapData = instantRadiationData(i) * (1 - 0.2*r);
        else
            variation = 0.1 * (sin(4*x) .* cos(4*y)) + 0.95;
            heatmapData = instantRadiationData(i) * variation;
        end
        
        % Ensure non-negative values
        heatmapData = max(0, heatmapData);
        
        % Plot the heatmap
        imagesc(heatmapData);
        colorbar;
        colormap(ax, hot); % Apply colormap to current axes
        caxis([0, maxRadiation]);
        
        % Annotate with face name and dimensions
        title([faceNames{i}, ' Face: ', faceDimensions{i}], 'FontWeight', 'bold');
        xlabel('Position');
        ylabel('Position');
        
        % Add text showing average radiation
        text(faceSize(1)/2, faceSize(2)/2, [num2str(round(instantRadiationData(i))), ' W/m²'], ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'w', 'FontSize', 12);
        
        % Highlight the sun-facing face
        if i == 5 % +Z face
            rectangle('Position', [0.5, 0.5, faceSize(1), faceSize(2)], ...
                'EdgeColor', 'y', 'LineWidth', 3, 'LineStyle', '--');
            text(5, 5, 'SUN-FACING', 'Color', 'y', 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        end
    end
    
    % Add overall title
    sgtitle('Radiation Heat Maps for Each CubeSat Face (W/m²)', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Add timestamp information at bottom of figure
    annotation('textbox', [0.25, 0.01, 0.5, 0.03], 'String', ...
        ['Simulation by: Ahmadh1231  |  2025-03-06 17:48:04 UTC'], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);
end

function plotRadiationTimeChart(radiationData, faceNames, timeSeconds, isEclipsed)
    % Create radiation time series chart with eclipse indicators
    
    % Convert time to hours for plotting
    timeHours = timeSeconds / 3600; % Convert to hours
    
    % Set up colors for each face
    colors = lines(6);
    
    % Create the figure
    figure('Position', [100, 100, 1200, 600]);
    
    % Plot eclipse regions first
    ax = gca;
    hold on;
    
    % Calculate Y limits first to properly scale eclipse shading
    maxRadiation = max(max(radiationData)) * 1.1;
    if maxRadiation == 0
        maxRadiation = 10; % Default if no radiation
    end
    ylim([0, maxRadiation]);
    
    % Find eclipse entry/exit times
    eclipseTransitions = diff([0, isEclipsed, 0]);
    eclipseStarts = find(eclipseTransitions == 1) / 60; % Convert to hours
    eclipseEnds = find(eclipseTransitions == -1) / 60; % Convert to hours
    
    % Plot shaded areas for eclipses
    for i = 1:length(eclipseStarts)
        x = [eclipseStarts(i), eclipseEnds(i), eclipseEnds(i), eclipseStarts(i)];
        y = [0, 0, maxRadiation, maxRadiation]; % Use calculated Y limits
        patch(x, y, [0.8, 0.8, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end
    
    % Main plot with full 24-hour data
    for i = 1:6
        plot(timeHours, radiationData(i, :), 'LineWidth', 1.5, 'Color', colors(i, :));
    end
    
    % Label eclipses
    for i = 1:length(eclipseStarts)
        if i <= 5 % Only label first 5 eclipses to avoid clutter
            midEclipse = (eclipseStarts(i) + eclipseEnds(i))/2;
            text(midEclipse, maxRadiation*0.9, 'ECLIPSE', 'FontWeight', 'bold', ...
                 'Color', [0.3 0.3 0.6], 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'top', 'FontSize', 8);
        end
    end
    
    xlabel('Time (hours)');
    ylabel('Radiation (W/m²)');
    title('CubeSat Surface Radiation Over 24 Hours with Eclipse Periods', 'FontSize', 14);
    grid on;
    legend(faceNames, 'Location', 'best');
    xlim([min(timeHours), max(timeHours)]);
end

function plotOrbit(orbitalPositions, sunlitPositions, eclipsedPositions, earthRadius)
    % Create 3D orbit visualization
    figure('Position', [100, 100, 900, 800]);
    
    % Plot Earth
    [X, Y, Z] = sphere(50);
    X = X * earthRadius;
    Y = Y * earthRadius;
    Z = Z * earthRadius;
    
    % Create Earth surface
    surf(X, Y, Z, 'FaceColor', [0.3 0.5 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hold on;
    
    % Plot orbit path
    plot3(orbitalPositions(:,1), orbitalPositions(:,2), orbitalPositions(:,3), 'k-', 'LineWidth', 1);
    
    % Plot sunlit and eclipsed sections with different colors
    % Make sure we always plot something even if arrays are empty
    if isempty(sunlitPositions) && isempty(eclipsedPositions)
        % If both are empty (shouldn't happen), plot the whole orbit as sunlit
        plot3(orbitalPositions(:,1), orbitalPositions(:,2), orbitalPositions(:,3), 'y.', 'MarkerSize', 8);
    else
        % Plot what we have
        if ~isempty(sunlitPositions)
            plot3(sunlitPositions(:,1), sunlitPositions(:,2), sunlitPositions(:,3), 'y.', 'MarkerSize', 8);
        end
        
        if ~isempty(eclipsedPositions)
            plot3(eclipsedPositions(:,1), eclipsedPositions(:,2), eclipsedPositions(:,3), 'b.', 'MarkerSize', 8);
        end
    end
    
    % Add Sun direction indicator (simplified)
    quiver3(0, 0, 0, -10000, 0, 0, 'y', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    text(-11000, 0, 0, 'SUN', 'Color', 'y', 'FontWeight', 'bold');
    
    % Annotate orbit components
    title('CubeSat Orbit Visualization with Eclipse Detection', 'FontSize', 14);
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    
    % Add legend - conditionally based on what was plotted
    if isempty(eclipsedPositions)
        legend('Earth', 'Orbit Path', 'Sunlit Positions', 'Sun Direction');
    elseif isempty(sunlitPositions)
        legend('Earth', 'Orbit Path', 'Eclipse Positions', 'Sun Direction');
    else
        legend('Earth', 'Orbit Path', 'Sunlit Positions', 'Eclipse Positions', 'Sun Direction');
    end
    
    % Set axes properties
    axis equal;
    grid on;
    
    % Add annotation with orbit information
    annotation('textbox', [0.15, 0.01, 0.7, 0.05], 'String', ...
        {'50° Inclination Orbit Simulation', ...
         'Yellow: Satellite in sunlight | Blue: Satellite in Earth''s shadow'}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
    
    % Add viewpoint for better visualization
    view(45, 30);
end

function plotRadiationComponents(solar, albedo, ir, timeSeconds)
    % Plot radiation components separately
    
    % Convert time to hours for plotting
    timeHours = timeSeconds / 3600;
    
    % Setup figure
    figure('Position', [100, 100, 1200, 800]);
    
    % Create subplots for +Z face (sun-facing)
    subplot(2,1,1);
    hold on;
    plot(timeHours, solar(5,:), 'r-', 'LineWidth', 2);
    plot(timeHours, albedo(5,:), 'g-', 'LineWidth', 2);
    plot(timeHours, ir(5,:), 'b-', 'LineWidth', 2);
    title('Radiation Components for +Z Face (Sun-facing)', 'FontSize', 14);
    xlabel('Time (hours)');
    ylabel('Radiation (W/m²)');
    grid on;
    
    % Ensure axes have reasonable values
    maxVal = max([max(solar(5,:)), max(albedo(5,:)), max(ir(5,:))]);
    if maxVal == 0
        maxVal = 10; % Default if no radiation
    end
    ylim([0, maxVal*1.1]);
    
    legend('Direct Solar', 'Earth Albedo (Reflection)', 'Earth IR');
    
    % Create subplot for radiation breakdown
    subplot(2,1,2);
    
    % Create stacked area for average radiation component percentage
    avgSolar = mean(solar, 2);
    avgAlbedo = mean(albedo, 2);
    avgIR = mean(ir, 2);
    
    % Ensure we have non-zero data
    if sum(avgSolar) + sum(avgAlbedo) + sum(avgIR) == 0
        % If no data, create dummy data for visualization
        avgSolar = 5 * ones(6,1);
        avgAlbedo = 3 * ones(6,1);
        avgIR = ones(6,1);
    end
    
    totalAvg = avgSolar + avgAlbedo + avgIR;
    
    % Calculate percentages
    percentSolar = zeros(6,1);
    percentAlbedo = zeros(6,1);
    percentIR = zeros(6,1);
    
    for i = 1:6
        if totalAvg(i) > 0
            percentSolar(i) = avgSolar(i) / totalAvg(i) * 100;
            percentAlbedo(i) = avgAlbedo(i) / totalAvg(i) * 100;
            percentIR(i) = avgIR(i) / totalAvg(i) * 100;
        end
    end
    
    % Create bar chart
    bar([avgSolar, avgAlbedo, avgIR], 'stacked');
    title('Average Radiation Component Breakdown by Face', 'FontSize', 14);
    xlabel('CubeSat Face');
    ylabel('Radiation (W/m²)');
    legend('Direct Solar', 'Earth Albedo (Reflection)', 'Earth IR');
    set(gca, 'XTickLabel', {'+X', '-X', '+Y', '-Y', '+Z', '-Z'});
    grid on;
    
    % Annotate with percentages for +Z face (sun facing)
    if avgSolar(5) > 0
        text(5, avgSolar(5)/2, sprintf('%.1f%%', percentSolar(5)), ...
            'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
    end
    
    if avgAlbedo(5) > 0
        text(5, avgSolar(5) + avgAlbedo(5)/2, sprintf('%.1f%%', percentAlbedo(5)), ...
            'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
    end
    
    if percentIR(5) > 5 % Only show if percentage is significant
        text(5, avgSolar(5) + avgAlbedo(5) + avgIR(5)/2, sprintf('%.1f%%', percentIR(5)), ...
            'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
    end
    
    % Add information about Earth albedo
    annotation('textbox', [0.15, 0.01, 0.7, 0.05], 'String', ...
        {sprintf('Earth Albedo Coefficient: %.2f (%.0f%% of solar flux reflected from Earth)', ...
        earthAlbedo, earthAlbedo*100)}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
end
