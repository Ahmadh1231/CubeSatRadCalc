% Run the CubeSat Radiation Simulation
% Created by: Ahmadh1231
% Date: 2025-03-06 15:29:01

fprintf('\n=== 6U CubeSat Radiation Simulation ===\n');
fprintf('Running simulation with 1-minute intervals over 24 hours...\n');

% Run the main simulation script
cubesat_radiation_simulation;

fprintf('\nSimulation completed.\n');
fprintf('Generated two output figures:\n');
fprintf('1. Radiation heatmap for each face\n');
fprintf('2. Radiation time chart over 24 hours\n');