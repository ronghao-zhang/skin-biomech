%% Script Information

% Script Introduction:
%     This script aims to create the stimulation trace use the downsample patterns. 

% Author:
%     Ronghao Zhang

% Dates:
%     Created: 2022-06-03
%     Last-Edited: 2022-08-23

% Methodology:
%    1. Create the stimulus trace matrix use the downsample pattern
%    2. Create stimulus and afferent response using the Stimulus()

%% Load data

% Use the following Command when using Bensmaia Lab PC
load('C:\Users\somlab\Desktop\skin-biomechanics\data\fabricScan\initialization_parameters.mat');
load('C:\Users\somlab\Desktop\skin-biomechanics\data\fabricScan\detrended_hmap.mat');

% remember to set up path for touchSim

%% Create the Trace for All Patterns

% Create an empty cell array to store the trace 
pattern_trace_cell = cell(2,num_patterns);

% Create the Trace Matrix
for n = 1:num_patterns
    
    % obtain the info from the previous cell array
    pattern_detrend = pattern_detrend_cell{2,n};
    width = pattern_detrend_cell{3,n};
    length = pattern_detrend_cell{4,n};
    adjust_height = pattern_detrend_cell{5,n};
    
    % create arrays that store the stripe x-y axis information
    axis_length = (pin_radius:touch_resolution:(length-1)*touch_resolution); % array that store the x-axis of the downsampled pattern
    axis_width = (pin_radius:touch_resolution:(width-1)*touch_resolution); % array that store the y-axis of the downsampled pattern height map
    
    % generate a padding to the height map
    pad_offset = padding*spat_freq; % number of pins need to be added before and after
    pattern_padding = ones(width,length+2*pad_offset)*adjust_height;
    
    for i = 1:length
        pattern_padding(:,(pad_offset+i)) = pattern_detrend(:,i);
    end
    
    % make new axis vectors for the padded strip
    pad_ax = (pin_radius:touch_resolution:(padding-pin_radius));
    axis_pad_length = [-fliplr(pad_ax), axis_length, pad_ax+axis_width(end)];
    
    % create the start point and the end point of scanning along the x axis
    [~, x_init] = min(abs(axis_pad_length + window_size/2)); 
    [~, x_term] = min(abs(axis_pad_length - (strip_length_array(1,n) + window_size/2)));

    % create the number of steps for the scanning
    n_steps = (x_term - x_init) / dt_idx; % Check if integer
    n_steps = floor(n_steps); 

    % create the y initial point (do not need to worry about it currently)
    y_init = 1; % [~, y_init] = min(abs(strip_axis_y + (window_size/2)));

    % set up index for the following for loop
    x_idx = x_init;
    y_idx = y_init;

    % set up the size for the trace that need to be input to the Stimulus()
    trace_input = zeros([n_steps,size(win_mesh,1)]); 

    % extrace the trace for each pin 
    for i = 1:n_steps

        trace = pattern_padding(y_idx : (y_idx+window_matrix_size-1), x_idx : (x_idx+window_matrix_size-1)); 
        trace(trace < 0) = 0;

        trace_input(i,:) = trace(:)';
        x_idx = x_idx + dt_idx;
    end

    % assign the values to the pattern_trace_cell
    pattern_trace_cell{1,n} = pattern_detrend_cell{1,n};
    pattern_trace_cell{2,n} = trace_input;
end

%% Create the Stimulus and Afferent Responses for all 14 patterns

% Create an empty cell array to store the stimulus
%   Row 1: name of the pattern
%   Row 2: stimulus 
%   Row 3: response
stim_resp_cell = cell(3,num_patterns);

% Create the Afferens Population
afferents = affpop_hand('D2d');

for n = 1:num_patterns
    
    trace = pattern_trace_cell{2,n};
    
    stim = Stimulus(trace, win_mesh, samp_freq, pin_radius); 
    resp = afferents.response(stim);
    
    stim_resp_cell{1,n} = pattern_trace_cell{1,n};
    stim_resp_cell{2,n} = stim;
    stim_resp_cell{3,n} = resp;
    
    disp(stim_resp_cell{1,n} + " is use to create stimulus and response successfully!")
    
end

%% Calculate the Indentation and Force

% Create a Cell to Store the Force by Time Plot
%   row1: the name of the pattern
%   row2: the time array for scanning
%   row3: the force array converted from indentation
force_time_cell = cell(3,num_patterns);

% Calculate the Conversion Parameter between Indentation and Force
force_per_indent_2mm = 0.05; % Unit: (Newton/mm indentation)/pin grid with 2mm diameter
previous_window_size = 2; % unit: mm, previous window area is 4mm^2.
previous_window_area = previous_window_size^2; % unit: mm^2
current_window_area = window_size^2;

num_pins = size(win_mesh,1); % the number of pins [!!!!!DOUBLE CHECK HERE!!!!!]

if num_pins ~= size(stim_resp_cell{2,1}.indentprofile,2)
    disp("ERROR: num_pins does not match indentprofile size!!!")
end

force_per_pin_indent = (force_per_indent_2mm/previous_window_area) * current_window_area/num_pins; % Unit: Newton/mm/pin

% This is the unit time for each movement
unit_time = dt_distance/scan_speed;

% Calculate the Force at each exact time point
for n = 1:num_patterns
    % obtain the indentation profile
    indent_pfl = stim_resp_cell{2,n}.indentprofile;
    
    % calculate the indentation by time and the applied force by time
    indent_by_time = sum(indent_pfl,2); % Unit: mm, each row is an total indentation, ordered by time sequence
    force_by_time = indent_by_time * force_per_pin_indent; % Unit: N, N/mm*mm, each row is the force applied on the texture by finger
    
    % calculate the time array
    time_array = zeros(size(force_by_time));
    
    % assign values to the force time cell
    force_time_cell{1,n} = stim_resp_cell{1,n};
    force_time_cell{2,n} = (0:unit_time:unit_time*(size(indent_pfl,1)-1))';
    force_time_cell{3,n} = force_by_time;
end

% calculate the average of the force applied for each pattern
force_avg = zeros(1,num_patterns);
for n = 1:num_patterns
    force_avg(1,n) = mean(force_time_cell{3,n},1);
end

%% Plot the Figures - Force v.s. Time
figure;
hold on
for n = 1:num_patterns/2
    subplot(2,7,n)
    plot(force_time_cell{2,n*2},force_time_cell{3,n*2})
    yline(force_avg(1,2*n),'color','#a61b29');
    ylim([0 0.4]);
    xlabel('time (sec)');
    ylabel('force (Newton)');
    title(force_time_cell{1,n*2});
end

for n = 1:num_patterns/2
    subplot(2,7,n+7)
    plot(force_time_cell{2,n*2-1},force_time_cell{3,n*2-1})
    yline(force_avg(1,2*n-1),'color','#a61b29');
    ylim([0 0.4]);
    xlabel('time (sec)');
    ylabel('force (Newton)');
    title(force_time_cell{1,n*2-1});
end
hold off

%% Plot The Response of Neurons 
figure
hold on
for n = 1:num_patterns/2
    subplot(2,7,n)
    hold on
    plot(stim_resp_cell{3,2*n});
    title("Response of " + pattern_detrend_cell{1,n*2});
    hold off
    
    subplot(2,7,n+7)
    hold on
    plot(stim_resp_cell{3,2*n-1})
    title("Response of " + pattern_detrend_cell{1,n*2-1});
    hold off
end
hold off

% Save the Stimulus and Response
save('C:\Users\somlab\Desktop\skin-biomechanics\data\fabricScan\afferent_response.mat', 'stim_resp_cell');

%% Calculate the Customized Pre-Indentation

force_expect = 0.25; % Unit: N, the expect force applied on each texture is 0.05N.
pre_indent_update = zeros(1,num_patterns);
for n = 1:num_patterns
    force_pre = force_avg(1,n); % current force applied for each pattern
    force_adj = force_expect - force_pre; % amount of force required to get 0.05N
    indent_adj = (force_adj/force_per_pin_indent) ... % Unit: N /(N/mm/pin), the indentation of just one pin to reach the amount of force
        /num_pins;
    pre_indent_update(1,n) = indent_adj + pre_indent;
end

%% Save the Stimulus and Response
save('C:\Users\somlab\Desktop\skin-biomechanics\data\fabricScan\pre_indentation.mat', 'pre_indent_update');
