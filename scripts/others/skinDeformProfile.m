%% Script Information

% Script Introduction:
%     This script aims to get the skin deformation information when the 
%     finger tip presses a certain square grating. 

% Author:
%     Ronghao Zhang

% Dates:
%     Created: 2022-06-15
%     Last-Edited: 2022-08-24

% Methodology:
%    1. Create Stimulus trace matrix for one unit time
%    2. Extract the indent Profile from there
%    3. Find the artifect location 

%% Initialization

% load the parameters
touch_resolution = 0.1;
spat_freq = 1/touch_resolution;
window_size = 10; 
pin_radius = touch_resolution/2;
samp_freq = 800;

% create a window matrix: size is 5mm*5mm (size need to times * spatial frequency)
window_matrix_size = window_size * spat_freq;
window_offset = window_matrix_size / 2; % half of 'pins' should be in negtive regions of xy-coordinate system
window_ax = linspace(-(window_size/2) + pin_radius,(window_size/2) - pin_radius, window_matrix_size);
[winx, winy] = meshgrid(window_ax,window_ax);
win_mesh = [winx(:) winy(:)];

%% Create a New Height map With Larger Size
% here use the 40th height map as an example
gap_width = 1.5; % gap_width (unit: mm) is defined as the width of the gap between gratings
grate_width = 0.5;  % grate_width (unit:mm) is defined as the width of the grating
grate_height = 1.0; 

gap_width_pix = gap_width*spat_freq;
grate_width_pix = grate_width*spat_freq;
unit_num = floor(window_matrix_size/(gap_width_pix + grate_width_pix));

col_height = zeros(1,window_matrix_size);

offset = floor(gap_width_pix/2);
idx_init = offset + 1; 
for i = 1:unit_num
    idx_term = idx_init + grate_width_pix - 1;
    col_height(1, idx_init:idx_term) = 1; 
    idx_init = idx_term + 1 + gap_width_pix;
end

col_height = col_height*grate_height;

height_map = ones(window_matrix_size);
height_map = height_map.*col_height; 

%% Obtain the IndentProfile from Stimulus

% create the trace for stimulation
trace =  zeros([1,size(win_mesh,1)]); 
trace(1,:) = height_map(:)';

% use the trace and window mesh to create stimulus
s = Stimulus(trace, win_mesh, samp_freq, pin_radius);

% create another empty matrix for skin indentation profile
skin_indent = zeros(size(height_map));

% assign the values from the indentprofile to the matrix
idx = 1; % the index for each row
d_col = window_matrix_size; % number of pixels in each height
for i = 1:d_col 
   skin_indent(i,:) = s.indentprofile(1,idx:idx+d_col-1); 
   idx = idx + d_col;
end

%% Figure: Plot the raw figure

% Create Scales for Imagesc Plots
pix_scale = linspace(touch_resolution,touch_resolution*window_matrix_size,window_matrix_size);
pix_scale = pix_scale - pin_radius;

% Comparing the Original HeightMap and Skin Indentprofile from the touchSim Package 
figure; 
hold on
subplot(1,2,1)
hold on
imagesc(pix_scale,pix_scale,height_map');
colormap hot;
title("The Height Map of the Non-Compliant Grating");
xlabel('Width (mm)')
ylabel('Length (mm)')
set(gca,'YDir','normal');
axis([0 10 0 10]);
hold off

subplot(1,2,2)
hold on
imagesc(pix_scale,pix_scale,skin_indent);
title("The Indentation Profile of the Grating Simulated Using touchSim");
xlabel('Width (mm)')
ylabel('Length (mm)')
set(gca,'YDir','normal');
axis([0 10 0 10]);
hold off
hold off

%% Reduce the Artifect

% get the loaction of gaps and indentations
col_height_nzero = find(col_height); % find the index of non-zero elements

col_height_gap = ones(size(col_height_nzero)); % empty matrix to store where "gaps exists"
for i = 1:size(col_height_nzero,2)-1
    if col_height_nzero(1,i+1) ==  col_height_nzero(1,i) + 1
        col_height_gap(1,i) = 0; % after 1, there should be gap
    end  
end
col_height_gap_nzero = find(col_height_gap);

% find the index where the gap starts and ends
gap_start = col_height_nzero(1, col_height_gap_nzero(1,1)) + 1;
gap_end = col_height_nzero(1, col_height_gap_nzero(1,1) + 1) - 1;

% find the index of the middle of gap --> where artifects is largest
mid_gap_idx = floor((gap_end - gap_start + 1)/2) + gap_start; 

temp = skin_indent(mid_gap_idx,:);

figure;
plot(1:window_matrix_size,temp)
