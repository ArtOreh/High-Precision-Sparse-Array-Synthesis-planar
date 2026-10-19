%% I. SYSTEM PARAMETERS AND DATA IMPORT

clear; clc; close all;

fc = 2.4e9;                     % Carrier frequency (Hz)
c = physconst("LightSpeed");    % Speed of light (m/s)
lambda = c / fc;                % Free-space wavelength (m)

% Maximum radiating aperture radius scaled by wavelength
Aperture_R_lambda = 5.0; 

% Additional structural boundary/clearance margin
board_clearance = 0.1 * lambda; 

% NOTE: The total physical radius of the structure is defined as:
% Total Radius = (Aperture_R_lambda * lambda) + Single_Patch_Radius + board_clearance

R_ground = Aperture_R_lambda*lambda; % Ground plane radius

filename = 'coordinates_2d_N=77_r=5.0lam0.85.txt'; % Input array configuration file

if exist(filename, 'file') == 2
    fileID = fopen(filename, 'r');
    fgetl(fileID); % Skip opening bracket or header line
    data = textscan(fileID, '[%f, %f],'); % Parse 2D coordinate pairs
    fclose(fileID);
    pts_normalized = [data{1}, data{2}]; % Format as Nx2 coordinate matrix
else
    error('Error: File %s not found. Ensure it is in the working directory.', filename);
end

N_elements = size(pts_normalized, 1);
fprintf('Successfully imported %d array elements.\n', N_elements);

pos_matrix = pts_normalized * R_ground; % Scale normalized points to physical dimensions

%% II. RADIATING ELEMENT SPECIFICATION AND DESIGN
% Circular microstrip patch design initialization
patch_template = patchMicrostripCircular;
patch_template.Substrate = dielectric('AIR'); 

% Generate baseline antenna dimensions at center frequency
auto_patch = design(patch_template, fc);

% Overwrite automated dimensions with optimized manual parameters
auto_patch.Radius = 35.2e-3;
auto_patch.FeedOffset = [10.7e-3, 0];
auto_patch.Height = 2.0e-3;

%% III. ARRAY ASSEMBLY VIA PCBSTACK
ant = pcbStack;

% Assign dielectric substrate thickness from the baseline design
ant.BoardThickness = auto_patch.Substrate.Thickness;

% Define the composite ground plane for Layer 1 (bottom layer)
% Considers core aperture, individual patch radius, and clearance margin
ground_layer = antenna.Circle('Radius', R_ground + auto_patch.Radius + board_clearance);
ant.BoardShape = ground_layer; 

% Construct Layer 2 (top layer) using the layout matrix and patch geometry
top_layer = antenna.Circle('Radius', auto_patch.Radius, 'Center', pos_matrix(1,:));
for i = 2:N_elements
    next_patch = antenna.Circle('Radius', auto_patch.Radius, 'Center', pos_matrix(i,:));
    top_layer = top_layer + next_patch; % Boolean union of top layer geometry
end

% Assign stack layers (top conductor, substrate, bottom ground plane)
ant.Layers = {top_layer, auto_patch.Substrate, ground_layer};

% Initialize feed location matrix
feed_locations = zeros(N_elements, 4);
for i = 1:N_elements
    f_x = pos_matrix(i,1) + auto_patch.FeedOffset(1);
    f_y = pos_matrix(i,2) + auto_patch.FeedOffset(2);
    
    % Define probe feeds transitioning vertically from Layer 3 (GND) to Layer 1 (Patch)
    feed_locations(i, :) = [f_x, f_y, 3, 1]; 
end
ant.FeedLocations = feed_locations;


% Generate and display the mesh layout (GrowthRate set to default 0.7 for preview)
meshData = mesh(ant, 'GrowthRate', 0.7); 

fprintf('  Number of triangles: %d\n', meshData.NumTriangles);
fprintf('  Number of Tetrahedra: %d\n', meshData.NumTetrahedra);
%% IV. VISUALIZATION AND ANALYSIS
% Render the physical layout of the array
figure('Name', 'Array Structure on Circular Ground Plane');
show(ant);
title(sprintf('Monolithic Ground Plane (R = %d mm) with Substrate Parameters', ground_layer.Radius*1e3));

% Plot the structural mesh
figure;
mesh(ant); 
title('Computational Mesh of the Patch Array');

% Evaluate computational resource requirements
memory_need = memoryEstimate(ant, fc);
fprintf('Estimated RAM requirements for EM simulation: %s\n\n', memory_need);
memory;

%% V. SYSTEM PARAMETER LOGGING AND VERIFICATION
fprintf('\n==================================================\n');
fprintf('                SYSTEM PARAMETERS                 \n');
fprintf('==================================================\n');

% Validate and log substrate characteristics
if ischar(auto_patch.Substrate.Name) || isstring(auto_patch.Substrate.Name)
    fprintf('Substrate Material Type:              %s\n', auto_patch.Substrate.Name);
else
    fprintf('Substrate Material Type:              Custom Dielectric\n');
end
fprintf('Relative Permittivity (epsR):         %.2f\n', auto_patch.Substrate.EpsilonR);
fprintf('Substrate Thickness:                  %.3f mm\n', auto_patch.Substrate.Thickness * 1e3);

fprintf('--------------------------------------------------\n');
fprintf('SINGLE RADIATING ELEMENT GEOMETRY:\n');
fprintf('Patch Radius:                         %.2f mm\n', auto_patch.Radius * 1e3);
fprintf('Patch Diameter:                       %.2f mm\n', (auto_patch.Radius * 2) * 1e3);
fprintf('Feed Point Offset (Relative to Patch Center):\n');
fprintf('  - X-axis Offset:                    %.2f mm\n', auto_patch.FeedOffset(1) * 1e3);
fprintf('  - Y-axis Offset:                    %.2f mm\n', auto_patch.FeedOffset(2) * 1e3);

fprintf('--------------------------------------------------\n');
fprintf('MONOLITHIC ARRAY GROUND PLANE PARAMETERS:\n');
fprintf('Ground Plane Geometry:                Circular\n');
fprintf('Composite Ground Radius (R_ground):   %.2f mm\n', ground_layer.Radius * 1e3);
fprintf('Composite Ground Diameter:            %.2f mm\n', (ground_layer.Radius * 2) * 1e3);

%% VI. FULL-WAVE MOM RADIATION PATTERN COMPUTATION AND DATA EXPORT

%return;

current_time = string(datetime('now', 'Format', 'dd.MM.yyyy HH:mm:ss'));
fprintf('\n[%s] Launching 3D Radiation Pattern computation. Computing MoM matrices...\n', current_time);

% Define high-resolution coordinate grid (1-degree resolution)
% Full upper-hemisphere coverage: Azimuth [-180, 180], Elevation [0, 90]
az_angles = -180:1:180;
el_angles = 0:1:90;

% Execute full-wave Method of Moments (MoM) field solver
[dir_values, az_grid, el_grid] = pattern(ant, fc, az_angles, el_angles);

% Export computed fields to a .mat file for seamless Python post-processing
save('mom_radiation_data_air.mat', 'az_grid', 'el_grid', 'dir_values');

fprintf('==================================================\n');
fprintf('     MOM SIMULATION DATA SUCCESSFULLY EXPORTED    \n');
fprintf('     Target File: mom_radiation_data_air.mat      \n');
fprintf('==================================================\n');

figure('Name', '3D Radiation Pattern View');
pattern(ant, fc); 

current_time = string(datetime('now', 'Format', 'dd.MM.yyyy HH:mm:ss'));
fprintf('Simulation completion timestamp: %s\n', current_time);

title('3D Radiation Pattern (MoM Solver Verified)');

%% VII Compute the return loss
% Compute the full multi-port scattering matrix at the center frequency
spar = sparameters(ant, fc);            

% Extract the raw complex S-parameter matrix (NxN complex numbers)
S_complex = spar.Parameters; 

% Calculate the Active S11 for each port (complex row sum for uniform excitation)
% This accounts for both geometric return loss and mutual coupling phases
active_S11_raw = sum(S_complex, 2); 

% Convert the active reflection coefficients to the logarithmic scale (dB)
active_S11_db = 20 * log10(abs(active_S11_raw));

% Display active reflection coefficients for all array elements in the console
fprintf('\n--------------------------------------------------\n');
fprintf('ACTIVE REFLECTION COEFFICIENTS (S11 ACTIVE) FOR ALL PORTS:\n');
fprintf('(Evaluated under simultaneous uniform array excitation)\n');
for idx = 1:N_elements
    fprintf('  - Port %3d: %.2f dB\n', idx, active_S11_db(idx));
end
fprintf('--------------------------------------------------\n');

% Calculate the physically correct mean active S11 using linear power scale
mean_active_S11_linear = mean(abs(active_S11_raw));
mean_active_S11_db = 20 * log10(mean_active_S11_linear);

fprintf('OVERALL ARRAY PERFORMANCE METRICS:\n');
fprintf('  - Mean Active S11:  %.2f dB\n', mean_active_S11_db);
fprintf('  - Worst-Case S11:   %.2f dB\n', max(active_S11_db));
fprintf('--------------------------------------------------\n');







