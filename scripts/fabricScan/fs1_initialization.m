%% Script Information

% Script Introduction:
%     This script specifies parameters that will be used to simulate afferent
%     response via touchSim Package when scanning over fabric surface using fingertips.

% Author:
%     Ronghao Zhang

% Dates:
%     Created: 2022-06-03
%     Last-Edited: 2022-08-23

% Methodology:
%    1. Define parameters, please see below for more details. 

%% Initialization

% Clear All On-Going Work
clc;
close all;

% Change the Directory to the Current Working Directory
cd('C:\Users\somlab\Desktop\skin-biomechanics\raw\fabric2021\fabric\')

% Set Path to the Old Patterns and Obtain the .csv Informations
old_patterns = dir(fullfile(pwd,'*.csv'));

% Parameters for Detrend the Desampled Height Map
window_size = 5;          % unit:mm, the width of the texture, the width & Length of the sensory area
touch_resolution = 0.1; %0.025 % unit:mm, the diamater of each 'pin' that contact the skin

if touch_resolution ~= 0.025
    fprintf('You are using a touch resolution of %.4f mm. It is OK but a normal resolution is 0.025 mm \n\n', touch_resolution);
end

% Parameters for Creating Stimulus and Response
scan_speed = 80; % unit: mm/s, the speed of scanning on texture (using finger)
samp_freq = 800; % unit: Hz, or number of pins being scanned per second.
pre_indent = 0.5; % unit: mm, the pre-indentation value added to all pins of the detrended height map
padding = 5; % unit: mm, the distance (delay) added to the pattern texture

% Parameters auto-generated
spat_freq = 1/touch_resolution; % unit: mm^(-1), Spatial frequency per mm
pin_radius = touch_resolution/2; % unit:mm, the radiux of each 'pin'
num_patterns = size(old_patterns,1); % the number of patterns in the directory
dt_distance = scan_speed / samp_freq; % unit: mm, the displacement after one time unit

dt_idx = dt_distance / touch_resolution; % Check if this is an intege

% Create Arrays to Store the Strip Resolutions
pat_resol_length = zeros(1,num_patterns); % store the x resolution of the strip
pat_resol_width = zeros(1,num_patterns);  % store the y resolution of the strip

for n = 1:num_patterns
    % assign the strip resolutions (both length and width)
    file_name = old_patterns(n).name;
    pat_resol_length(n) = dlmread(file_name,',',[3,1,3,1])/1000; % unit: um/1000 -> mm per pixel, pixel size in length direction
    pat_resol_width(n) = dlmread(file_name,',',[4,1,4,1])/1000;  % unit: um/1000 -> mm per pixel, pixel size in width direction
end

% Store the Names of the Patterns
pattern_order = ["Chiffon (Gel)"         "Chiffon (No Gel)"         "Denim (Gel)"        "Denim (No Gel)"       ...
                 "Velveteen (Gel)"       "Velveteen (No Gel)"       "Microsuede (Gel)"   "Microsuede (No Gel)"  ...
                 "Satin (Gel)"           "Satin (No Gel)"           "Snowflake (Gel)"    "Snowflake (No Gel)"   ...
                 "Wool Gabardine (Gel)"  "Wool Gabardine (No Gel)"                                              ...
                 ];

% Create Window Size
window_matrix_size = window_size * (1/touch_resolution); % 1/res is the number of 'pins' per mm distance
window_offset = window_matrix_size / 2; % half of 'pins' should be in negtive regions of xy-coordinate system
window_ax = linspace(-(window_size/2) + pin_radius,(window_size/2) - pin_radius, window_matrix_size);
[winx, winy] = meshgrid(window_ax,window_ax);
win_mesh = [winx(:) winy(:)];

%% Save All Works in This Script

% Delete Unnecessary Variables 
clear n 

% Output the Data
save('C:\Users\somlab\Desktop\skin-biomechanics\data\fabricScan\initialization_parameters.mat');
    