%% 6U CubeSat Radiation Simulation
% This script simulates a 6U CubeSat with sun pointing on its largest side (+Z)
% and generates radiation heatmaps and time charts
% Created by: Ahmadh1231
% Date: 2025-03-06 15:29:01

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
earthAlbedo = 0.3; % Earth's albedo coefficient
earthIR = 237; % W/m^2 (Earth's IR radiation)

% Orbit parameters
altitude = 500; % km
earthRadius = 6371; % km
orbitRadius = earthRadius + altitude;

% Material properties
absorptivity = 0.8; % Solar absorptivity
emissivity = 0.8; % Infrared emissivity

% Simulation parameters
simulationHours = 24;  % 24-hour simulation
minutesPerHour = 60;
timeStepMinutes = 1;   % 1-minute intervals
numTimeSteps = uint32(simulationHours * minutesPerHour / timeStepMinutes); % Integer type conversion
orbitPeriod = 2 * pi * sqrt(orbitRadius^3 / (3.986e14)); % seconds

% Faces 
numFaces = 6;

%% Face data
faceNames = {'+X (2U)', '-X (2U)', '+Y (3U)', '-Y (3U)', '+Z (6U)', '-Z (6U)'};
faceDimensions = {'10cm × 20cm', '10cm × 20cm', '10cm × 30cm', '10cm × 30cm', '20cm × 30cm', '20cm × 30cm'};

%% Create time and results arrays with guaranteed integer size
timeSeconds = zeros(1, numTimeSteps);
radiationData = zeros(numFaces, numTimeSteps);
earthViewFactors = zeros(numFaces, numTimeSteps);

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

%% Main simulation loop
for t = 1:numTimeSteps
    % Current simulation time in seconds
    currentTime = timeSeconds(t);
    
    % Orbital position (angle in radians)
    theta = 2 * pi * mod(currentTime / orbitPeriod, 1);
    
    % CubeSat position in orbit (assuming circular orbit)
    position = orbitRadius * [cos(theta), sin(theta), 0];
    
    % Sun position (assuming sun is in +Z direction)
    % Sun position changes over 24 hours to simulate Earth's rotation
    sunOffset = 2 * pi * (currentTime / (24 * 3600)); 
    sunPosition = [sin(sunOffset), 0, cos(sunOffset)] * 1e8; 
    
    % Calculate sun direction
    sunDirection = sunPosition - position;
    sunDirection = sunDirection / norm(sunDirection);
    
    % Required rotation to point +Z towards sun
    defaultZAxis = [0, 0, 1];
    rotationAxis = cross(defaultZAxis, sunDirection);
    
    % Create rotation matrix
    if norm(rotationAxis) < 1e-10
        rotationMatrix = eye(3);
    else
        rotationAxis = rotationAxis / norm(rotationAxis);
        rotationAngle = acos(dot(defaultZAxis, sunDirection));
        
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
    faceNormals = zeros(numFaces, 3);
    for i = 1:numFaces
        face = faces(i, :);
        v1 = rotatedVertices(face(2), :) - rotatedVertices(face(1), :);
        v2 = rotatedVertices(face(3), :) - rotatedVertices(face(2), :);
        normal = cross(v1, v2);
        faceNormals(i, :) = normal / norm(normal);
    end
    
    % Calculate nadir direction
    nadir = -position / norm(position);
    
    % Earth rho angle
    rho = asin(earthRadius / (earthRadius + altitude));
    
    % Calculate radiation for each face
    for faceIdx = 1:numFaces
        % Calculate radiation on each face
        
        % 1. Earth view factor
        angle = acos(dot(faceNormals(faceIdx,:), nadir));
        if angle < pi/2 && angle < rho
            % Face can see Earth
            earthViewFactors(faceIdx, t) = cos(angle) * (1 - cos(rho)) / 2;
        else
            % Face cannot see Earth
            earthViewFactors(faceIdx, t) = 0;
        end
        
        % 2. Solar radiation
        cosIncidence = dot(faceNormals(faceIdx,:), sunDirection);
        if cosIncidence > 0
            solarRadiation = solarFlux * cosIncidence;
        else
            solarRadiation = 0; % Face doesn't see the sun
        end
        
        % 3. Earth albedo radiation
        albedoRadiation = solarFlux * earthAlbedo * earthViewFactors(faceIdx, t);
        
        % 4. Earth IR radiation
        irRadiation = earthIR * earthViewFactors(faceIdx, t);
        
        % 5. Total radiation
        radiationData(faceIdx, t) = solarRadiation + albedoRadiation + irRadiation;
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

%% 2. Radiation Time Chart
figure(2);
plotRadiationTimeChart(radiationData, faceNames, timeSeconds);

%% Helper Functions for Visualization

function plotRadiationHeatmaps(radiationData, instantRadiationData, faceNames, faceDimensions)
    % Create radiation heatmaps
    
    % Find global max radiation for consistent color scale
    maxRadiation = max(radiationData(:)) * 1.1;
    
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
            rectangle('Position', [0, 0, faceSize(1), faceSize(2)], ...
                'EdgeColor', 'y', 'LineWidth', 3, 'LineStyle', '--');
            text(5, 5, 'SUN-FACING', 'Color', 'y', 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        end
    end
    
    % Add overall title
    sgtitle('Radiation Heat Maps for Each CubeSat Face (W/m²)', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Add timestamp information at bottom of figure
    annotation('textbox', [0.25, 0.01, 0.5, 0.03], 'String', ...
        ['Simulation by: Ahmadh1231  |  ', '2025-03-06 15:29:01 UTC'], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);
end

function plotRadiationTimeChart(radiationData, faceNames, timeSeconds)
    % Create radiation time series chart
    
    % Convert time to hours for plotting
    timeHours = timeSeconds / 3600; % Convert to hours
    
    % Set up colors for each face
    colors = lines(6);
    
    % Create the figure
    figure('Position', [100, 100, 1200, 600]);
    
    % Main plot with full 24-hour data
    hold on;
    
    for i = 1:6
        plot(timeHours, radiationData(i, :), 'LineWidth', 1.5, 'Color', colors(i, :));
    end
    
    xlabel('Time (hours)');
    ylabel('Radiation (W/m²)');
    title('CubeSat Surface Radiation Over 24 Hours (1-minute intervals)', 'FontSize', 14);
    grid on;
    legend(faceNames, 'Location', 'best');
    xlim([min(timeHours), max(timeHours)]);

end