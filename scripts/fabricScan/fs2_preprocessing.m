%% Script Information

% Script Introduction:
%     This script perform pre-processing of the height map of 

% Author:
%     Ronghao Zhang

% Dates:
%     Created: 2022-06-03
%     Last-Edited: 2022-08-23

% Methodology:
%    1. Load the original height map
%    2. Downsample and detrend the height map
%    3. Plot the downsampled and detrended height map

%% Load the Patterns

% load the initialization file that contains user-defined parameters 
load('C:\Users\somlab\Desktop\skin-biomechanics\data\fabricScan\initialization_parameters.mat');

% load the customized pre-indentation
load('C:\Users\somlab\Desktop\skin-biomechanics\data\fabricScan\pre_indentation.mat');

% change the direcotry to the pattern directory
cd('C:\Users\somlab\Desktop\skin-biomechanics\raw\fabric2021\fabric\');

%% Downsample & Store the Height Map
 
% Creat Cell Array for Storing the downsampled Height Map 
pattern_downsample_cell = cell(4,num_patterns); % row1:filename, row2:downsample pattern, row3:length, row4:width

% Creat An Array for Storing the Actual Length of Each
strip_length_array = zeros(1,num_patterns);

% Calculate the Downsampled Height Map (this is version 01 Downsample Method)
for n = 1:num_patterns
    
    % define the current file name
    file_name = old_patterns(n).name;
    
    % load the raw pattern and corresponding size
    input_table = readtable(file_name,'PreserveVariableNames',true);
    raw_width_pix = size(input_table,1); % unit: num of pixels
    raw_length_pix = size(input_table,2) - 1; % unit: num of pixels, -1: the last row is empty
    pattern_raw = table2array(input_table(:,2:raw_length_pix)); % convert the input_table to matrix and delete the last column
    
    % calculate the number of pixel needed to create the strip
    strip_width_pixels = window_size/pat_resol_width(n);
    
    % create an offset for centralize the scanning
    raw_width_offset = round(raw_width_pix/2);
    width_init = raw_width_offset - round(strip_width_pixels/2) + 1; % index of the row that the strip extraction starts
    width_term = width_init + ceil(strip_width_pixels) - 1; % index of the row that the strip extraction ends -> given a stripe with size strip_width
    
    % obtain the down-sized sample
    pattern_dsz = pattern_raw(width_init:width_term,:);
    pat_dsz_width = size(pattern_dsz,1);
    pat_dsz_length = size(pattern_dsz,2);
    
    % calculate the stripe size
    strip_width = window_size; % unit: mm, actual width of the stripe
    strip_length = pat_resol_length(n) * pat_dsz_length; % unit: mm, actual length of the stripe
    strip_size = floor([strip_width * spat_freq, strip_length * spat_freq]); % use number of pins to represent length and width
    
    % Downsample Pattern I: set-up a downsampled stripe
    pattern_downsample = zeros([strip_size(1) strip_size(2)]); % create an empty downsampled matrix
    pix_pin_length = touch_resolution/pat_resol_length(n); % unit: mm per pin /(mm per pix)
    pix_pin_width = touch_resolution/pat_resol_width(n); % unit: mm per pin /(mm per pix)
    
    for i = 1:strip_size(1)
        % define the start row and end row of downsampling method
        i_init = floor((i-1)*pix_pin_width + 1);
        i_init_unrnd = (i-1)*pix_pin_width + 1;
        i_term = ceil(i*pix_pin_width);
        i_term_unrnd = i*pix_pin_width;
        
        % create the weight array for width
        [current_size_i, current_weight_i] = fitWeightArray(pix_pin_width,i_init,i_term,i,true);
        
        for j = 1:strip_size(2)
            
            % define the start col and end col of downsampling method
            j_init = floor((j-1)*pix_pin_length + 1);
            j_init_unrnd = (j-1)*pix_pin_length + 1;
            j_term = ceil(j*pix_pin_length);
            j_term_unrnd = j*pix_pin_length;
            
            % create the weight array for length
            [current_size_j, current_weight_j] = fitWeightArray(pix_pin_length,j_init,j_term,j,false);
            
            % create the pixel grid that need to be downsampled
            current_pix_grid = zeros(current_size_i,current_size_j);
            current_pix_grid(1:current_size_i,1:current_size_j) = pattern_dsz(i_init:i_term, j_init:j_term);
            
            % create the weight matrix 
            current_weight = ones(current_size_i,current_size_j);
            current_weight = current_weight.*current_weight_i;
            current_weight = current_weight.*current_weight_j;
            
            % calculate the weighted average of the current pixel grid
            current_pix_grid = current_pix_grid.*current_weight;
            pattern_downsample(i,j) = sum(current_pix_grid,'all')/sum(current_weight,'all');
            
        end
    end
    
    % touchSim will should take the height with units of mm, but the height map gives unit in um
    pattern_downsample = pattern_downsample/1000; % um/1000 = mm
    
    % store the file_name and the downsampled pattern to the cell and their size
    pattern_downsample_cell{1,n} = pattern_order(1,n); % name
    pattern_downsample_cell{2,n} = pattern_downsample; % downsampled pattern
    pattern_downsample_cell{3,n} = strip_size(1); %size(pattern_downsample,1);    % downsampled pattern width
    pattern_downsample_cell{4,n} = strip_size(2); %size(pattern_downsample,2);    % downsampled pattern length
    
    strip_length_array(1,n) = strip_size(2)*touch_resolution; % store the length of each array
    
    % display the progress
    disp(pattern_downsample_cell{1,n} + " is loaded and downsampled successfully!")
    
end

%% Detrending the Downsampled Height Maps

% create a cell array to store the the detrend information
%   row 1: name of the pattern
%   row 2: detrended pattern
%   row 3: size @ width side
%   row 4: size @ length side
%   row 5: store the adjust height
pattern_detrend_cell = cell(5, num_patterns);

% create another cell array to store the info for detrend plots
%   row 1: name of the pattern
%   row 2: length avg before detrend
%   row 3: length avg after detrend
%   row 4: width avg before detrend
%   row 5: width avg after detrend
detrend_plot_cell = cell(5, num_patterns);

for n = 1:num_patterns
    
    % obtain the width and length of the current downsampled height map
    pattern_downsample = pattern_downsample_cell{2,n};
    dt_width = pattern_downsample_cell{3,n};
    dt_length = pattern_downsample_cell{4,n};
    
    % calculate the length avg before detrending, after detrending, and the value need to be subtracted
    [length_avg_bef, length_dev] = fitModelOnAvgDt(pattern_downsample,dt_width,true);
    
    % peform detrending on each row (width side)
    pattern_detrend = zeros(size(pattern_downsample));
    length_avg_aft = zeros(size(length_avg_bef));
    for wid = 1:dt_width
        pattern_detrend(wid,:) = pattern_downsample(wid,:) - length_dev(wid);
        length_avg_aft(1,wid) = mean(pattern_detrend(wid,:));
    end
    
    % calculate the width avg before detrending, after detrending, and the value need to be subtracted
    [width_avg_bef, width_dev] = fitModelOnAvgDt(pattern_detrend,dt_length,false);
    
    width_avg_aft = zeros(size(width_avg_bef));
    for len = 1:dt_length
        pattern_detrend(:,len) = pattern_detrend(:,len) - width_dev(len);
        width_avg_aft(1,len) = mean(pattern_detrend(:,len));
    end
    
    % add the pre-indentation
    curr_min_height = prctile(pattern_detrend,5,"All"); % unit: um
    %%%%%%%%%%%%%%%%Change for customerize Indentation%%%%%%%%%%%%%%%%
    % adjust_height = abs(pre_indent - curr_min_height);
    adjust_height = abs(pre_indent_update(1,n) - curr_min_height);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pattern_detrend = pattern_detrend + adjust_height;
    
    % store the values for subsequent analysis
    pattern_detrend_cell{1,n} = pattern_order(1,n); % name
    pattern_detrend_cell{2,n} = pattern_detrend; % downsampled pattern
    pattern_detrend_cell{3,n} = pattern_downsample_cell{3,n};   % detrend pattern width remains the same
    pattern_detrend_cell{4,n} = pattern_downsample_cell{4,n};   % detrend pattern length remains the same
    pattern_detrend_cell{5,n} = adjust_height;
    
    % store the values for plotting the detrending process
    detrend_plot_cell{1,n} = pattern_order(1,n); % name
    detrend_plot_cell{2,n} = [1:dt_width ; length_avg_bef];
    detrend_plot_cell{3,n} = [1:dt_width ; length_avg_aft];
    detrend_plot_cell{4,n} = [1:dt_length ; width_avg_bef];
    detrend_plot_cell{5,n} = [1:dt_length ; width_avg_aft];
    
end

%% save the downsampled pattern and the length of each strip
save('C:\Users\somlab\Desktop\skin-biomechanics\data\fabricScan\detrended_hmap.mat', ... 
     'pattern_downsample_cell','strip_length_array','pattern_detrend_cell','detrend_plot_cell','strip_length_array');
 
%% Plot the Downsampled Height Map
figure 
hold on
for n = 1:num_patterns/2
    
    % use a 4 by 4 subplot system for plotting
    subplot(2,7,n)
    hold on
    % plot the texture pattern
    fabric = pattern_downsample_cell{2,n*2}';
    fabric_x = 1:touch_resolution:touch_resolution*size(fabric,2); 
    fabric_y = 1:touch_resolution:touch_resolution*size(fabric,1); 
    imagesc(fabric_x,fabric_y,fabric);
    
    % plot the texture size in mm
    xlim([1 max(fabric_x)]);
    ylim([1 max(fabric_y)]);
    xlabel("width(mm)");
    ylabel("length(mm)");
    
    % add the title and inverse the y-axis
    title("D.S. " + pattern_downsample_cell{1,n*2});
    
    % set up a colormap and present the color bar
    colormap hot;
    colorbar;
    
    % draw the color bar max and min
    cmax = max(max(pattern_downsample_cell{2,n*2}));
    cmin = min(min(pattern_downsample_cell{2,n*2}));
    caxis([cmin cmax]); 
    hold off
    
    subplot(2,7,n+7)
    hold on
    % plot the texture pattern
    fabric = pattern_downsample_cell{2,n*2-1}';
    fabric_x = 1:touch_resolution:touch_resolution*size(fabric,2); 
    fabric_y = 1:touch_resolution:touch_resolution*size(fabric,1); 
    imagesc(fabric_x,fabric_y,fabric);
    
    % plot the texture size in mm
    xlim([1 max(fabric_x)]);
    ylim([1 max(fabric_y)]);
    xlabel("width(mm)");
    ylabel("length(mm)");
    
    % add the title and inverse the y-axis
    title("D.S. " + pattern_downsample_cell{1,n*2-1});
    
    % set up a colormap and present the color bar
    colorbar;
    caxis([cmin cmax]); 
    hold off
end
hold off

%% Plot the Detrended Height Map

figure 
hold on
for n = 1:num_patterns/2
    
    % use a 2 by 7 subplot system for plotting
    subplot(2,7,n)
    hold on
    
    % plot the texture pattern
    fabric = pattern_detrend_cell{2,n*2}';
    fabric_x = 1:touch_resolution:touch_resolution*size(fabric,2); 
    fabric_y = 1:touch_resolution:touch_resolution*size(fabric,1); 
    imagesc(fabric_x,fabric_y,fabric);
    
    % plot the texture size in mm
    xlim([1 max(fabric_x)]);
    ylim([1 max(fabric_y)]);
    xlabel("width(mm)");
    ylabel("length(mm)");
    
    % add the title and inverse the y-axis
    title("D.T. " + pattern_detrend_cell{1,n*2});
    set(gca,'YDir','normal');
    
    % set up a colormap and present the color bar
    colormap hot;
    colorbar;
    
    % draw the color bar max and min
    cmax = max(max(pattern_detrend_cell{2,n*2}));
    cmin = min(min(pattern_detrend_cell{2,n*2}));
    caxis([cmin cmax]); 
    
    hold off
    
    
    subplot(2,7,n+7)
    hold on
    % plot the texture pattern
    fabric = pattern_detrend_cell{2,n*2-1}';
    fabric_x = 1:touch_resolution:touch_resolution*size(fabric,2); 
    fabric_y = 1:touch_resolution:touch_resolution*size(fabric,1); 
    imagesc(fabric_x,fabric_y,fabric);
    
    % plot the texture size in mm
    xlim([1 max(fabric_x)]);
    ylim([1 max(fabric_y)]);
    xlabel("width(mm)");
    ylabel("length(mm)");
    
    % add the title and inverse the y-axis
    title("D.T. " + pattern_detrend_cell{1,n*2-1});
    set(gca,'YDir','normal');
    
    % set up a colormap and present the color bar
    colormap hot;
    colorbar;
    caxis([cmin cmax]); 
    hold off
end
hold off


%% Function: Fit the Weights to Each Matrix
function [curr_size,current_weight] = fitWeightArray(pix_per_pin,init,term,idx,isWidth)
    
    % obtain the current size of the weight array
    curr_size = term - init + 1;
    
    % create the array (h array for length side, v array for width side)
    if isWidth
        current_weight = ones(curr_size,1);
    else 
        current_weight = ones(1,curr_size);
    end
    
    % obtain the start width
    weight_str = ceil((idx-1)*pix_per_pin) - (idx-1)*pix_per_pin;
    
    if weight_str == 0
        current_weight(1,1) = 1;
    else
        current_weight(1,1) = weight_str;
    end
    
    % obtain the end width
    weight_end = idx*pix_per_pin - floor(idx*pix_per_pin);
    
    if weight_end == 0
        current_weight(end) = 1;
    else
        current_weight(end) = weight_end;
    end
end

%% Function: Fit Model on Row or Column Average for Detrended Purposes
function [avg_bef, deviation] = fitModelOnAvgDt(pattern,oppo_size,is_length) % if is length, oppo_size should be width size
    
    % create an empty matrix to store the average value of each rol or col before detrending
    avg_bef = zeros(1,oppo_size); 
    
    % calculate the average of row or col based on isRow value
    if is_length 
        for len = 1:oppo_size
            avg_bef(len) = mean(pattern(len,:));
        end
    else 
        for col = 1:oppo_size
            avg_bef(col) = mean(pattern(:,col));
        end
    end
    
    % fit two models for the calculated average
    dep_input = (1:oppo_size)';
    indep_input = avg_bef';
    [fit_poly1,gof_poly1] = fit(dep_input,indep_input,"poly1");
    [fit_poly2,gof_poly2] = fit(dep_input,indep_input,"poly2");
   
    % select the model based on adjusted r^2 value
    if gof_poly1.adjrsquare > gof_poly2.adjrsquare
        deviation = polyval(coeffvalues(fit_poly1),1:oppo_size); 
    else
        deviation = polyval(coeffvalues(fit_poly2),1:oppo_size);
    end
end