%% Script Information

% Script Introduction:
%     This script aims to load, downsample, detrend the Gelsight of the non
%     compliant ojects. 

% Author:
%     Ronghao Zhang

% Dates:
%     Created: 2022-07-13
%     Last-Edited: 2022-08-05

% Methodology:
%    1. Load the SquareGrating m1 Gel Sights and Corresponding None-Gel Version
%    2. Calibrate, Detrending, Downsampling the GelSight
%    3. Detrending the Ground Truth
%    4. Cross Correlation 

% Check before Run
% 1. Replace all m? with expected gel index
% 2. directory of input and directory of output

%% Load the Gel and Non-Gel Data

order = ["Denim", "Microsuede", "SandPaper"]; 

% Slope and Intercept for Calibration
load('C:\Users\somlab\Desktop\skin-biomechanics\data\freqAnalysis\calibration\cal_linear_models.mat');
m1_slope = Slope_intcpt{1,1}(2);     %%%% NEED TO CHANGE HERE Slope_intcpt{1,1}=m_one {1,2}=m_three {1,3}=m_four {1,4}=m_five 
m1_intercept = Slope_intcpt{1,1}(1); %%%% NEED TO CHANGE HERE 

disp("Loading of Regression Model for Calibration is Completed!"); 

% Change the Directory to the Data Directory
cd('C:\Users\somlab\Desktop\skin-biomechanics\raw\gelsight2022\squareGratings\m1\')

% Obtain m1 GelSight and None GelSight Info
m1_gelsight = dir(fullfile(pwd,'.\m1*.csv'));
m1_grdtruth = dir(fullfile(pwd,'.\None*.csv'));

if size(m1_gelsight) == size(m1_gelsight) 
    m1_pair_num = size(m1_gelsight,1);
end

% Empty Data Structure to Store Info
m1_gs_res = zeros(2,m1_pair_num); % gelsight heightmap resolution x,y
m1_gt_res = zeros(2,m1_pair_num); % ground truth heightmap resolution x,y
m1_gs_hmap = cell(1,m1_pair_num); % gelsight heightmap
m1_gt_hmap = cell(1,m1_pair_num); % ground truth heightmap
m1_gs_size = zeros(2,m1_pair_num);  % gelsight heighmap size
m1_gt_size = zeros(2,m1_pair_num);  % ground truth heightmap size

for n = 1:m1_pair_num
   
   % --- for the gelsight ---
   % file names and folder path
   f_name = m1_gelsight(n).name; 
   f_path = m1_gelsight(n).folder; 
   
   % import the gel height map
   [res_x, res_y, hmap] = ImportHeightMap(f_name, f_path); 
    
   % export resolution, heightmap, size
   m1_gs_res(1,n) = res_x;
   m1_gs_res(2,n) = res_y;
   m1_gs_hmap{1,n} = hmap;
   m1_gs_size(1,n) = size(hmap,2);
   m1_gs_size(2,n) = size(hmap,1);
   
   % --- for ground truth ---
   % file names and folder path
   f_name = m1_grdtruth(n).name; 
   f_path = m1_grdtruth(n).folder; 
   
   % import the gel height map
   [res_x, res_y, hmap] = ImportHeightMap(f_name, f_path); 
    
   % export resolution, heightmap, size
   m1_gt_res(1,n) = res_x;
   m1_gt_res(2,n) = res_y;
   m1_gt_hmap{1,n} = hmap;
   m1_gt_size(1,n) = size(hmap,2);
   m1_gt_size(2,n) = size(hmap,1);
   
end

clear f_name f_path hmap res_x res_y
disp("Loading of Height Maps is Completed!");

%% Detrending the Ground Truth Image

% Empty cell array to store the detrended ground-truth height map
m1_gt_hmap_dt = cell(1,m1_pair_num); 

for n = 1:m1_pair_num
        
    % Detrending
    hmap_dt = Detrending(m1_gt_hmap{1,n});
    
    % Add Indentation Using 5% 
    prctile5_hgt = prctile(hmap_dt,5,'All');
    hmap_dt = hmap_dt - prctile5_hgt; 
    hmap_dt(hmap_dt<0) = 0;
    
    % export the hmap
    m1_gt_hmap_dt{1,n} = hmap_dt'; 
    
end

clear hmap_dt prctile5_hgt
disp("Detrending of the Ground Truth Height Map is Completed!");
%% Downsampling, Detrending, And Calibrating GelSight Height Map

% Empty Cell Array to Store the Detrended Ground-truth Height Map
m1_gs_hmap_dt = cell(1,m1_pair_num);   

for n = 1:m1_pair_num
    
    % obtain the height map and the dimension
    hmap = m1_gs_hmap{1,n}; 
    hmap_x = m1_gs_size(1,n);
    hmap_y = m1_gs_size(2,n);
    raw_x = m1_gt_size(1,n);
    raw_y = m1_gt_size(2,n);
    
    % calculate the downsampled height map
    hmap_ds = DownSampling(hmap,hmap_x,hmap_y,raw_x, raw_y);
    
    % calculate the detrended height map
    hmap_dt = Detrending(hmap_ds); 
    
    % Add Indentation Using 5% 
    prctile5_hgt = prctile(hmap_dt,5,'All');
    hmap_dt = hmap_dt - prctile5_hgt; 
    hmap_dt(hmap_dt<0) = 0;
    
    % Calibrate the Gel Sight Height Map
    hmap_dt = (hmap_dt - m1_intercept)/m1_slope;
    
    % Export and Store
    m1_gs_hmap_dt{1,n} = hmap_dt'; 
    
end

clear hmap hmap_x hmap_y raw_x raw_y hmap_ds hmap_dt prctile5_hgt hmap_dt
disp("Detrending and Calibration of the Gel Sight Heigh Map is Completed!");

%% Cross Correlation to find the Angle Rotation and Scaling Coefficient

% define the input parameter
min_degree = -3.5;
max_degree = 3.5;
deg_step_size = 0.1;
min_scale =  0.87;
max_scale = 0.95;
scale_step_size = 0.01;
chop_size = 500;

% find the best fit angle and scaling factor.
m1_cc_param = zeros(m1_pair_num,5);
m1_cc_output = cell(m1_pair_num,1);

for n = 1:m1_pair_num
    gel_hmap = m1_gs_hmap_dt{1,n};
    gt_hmap = m1_gt_hmap_dt{1,n};
    [max_cc, degree, scale_coef, offset, cc_matrix] = CalcCorrlation(gel_hmap,gt_hmap, ...
                                                                     min_degree, max_degree,deg_step_size, ...
                                                                     min_scale, max_scale, scale_step_size, ...
                                                                     chop_size);
    m1_cc_param(n,:) = [max_cc, degree, scale_coef, offset];
    m1_cc_output{n,1} = cc_matrix;
end

disp("Cross Correlation is Completed!");
%%
% adjust the height map for matching
m1_hmap_match = cell(m1_pair_num,2);
for n = 1:m1_pair_num
    gel_hmap = m1_gs_hmap_dt{1,n};
    gt_hmap = m1_gt_hmap_dt{1,n};
    
    degree = m1_cc_param(n,2);
    scale_coef = m1_cc_param(n,3);
    offset = m1_cc_param(n,4:5);
    
    gel_hmap_r = imrotate(gel_hmap,degree,'bilinear');
    gel_hmap_c = gel_hmap_r(chop_size+1:end-chop_size,chop_size+1:end-chop_size);
    gel_hmap_s = imresize(gel_hmap_c,scale_coef);
    adj_gel_hmap = gel_hmap_s;
    
    adj_gel_hmap_x = size(adj_gel_hmap,2);
    adj_gel_hmap_y = size(adj_gel_hmap,1);
    
    adj_gt_hmap = gt_hmap(offset(2)+1:offset(2)+adj_gel_hmap_y, offset(1)+1:offset(1)+adj_gel_hmap_x);
    
    m1_hmap_match{n,1} = adj_gt_hmap;
    m1_hmap_match{n,2} = adj_gel_hmap;
end

%% Save The Matched Height Map After Pre-Processing
save('C:\Users\somlab\Desktop\nonCompliantTouchSim_rz_2022\data\fabric\m1_fabric_preProcessing.mat'  , ...
     'order', 'm1_hmap_match','m1_cc_param','m1_cc_output', 'min_degree', 'max_degree', 'deg_step_size', ...
     'min_scale', 'max_scale', 'scale_step_size', 'chop_size','m1_gt_hmap_dt','m1_gs_hmap_dt');

%% Function: to import the resolution and Height Map
function [res_x, res_y, matrix] = ImportHeightMap(file_name,file_path)
    % set up path
    cd(file_path);
    % get the resolution
    res_x = dlmread(file_name,',',[3,1,3,1])/1000; % unit: mm per pixel, pixel size in length direction
    res_y = dlmread(file_name,',',[4,1,4,1])/1000;  % unit: mm per pixel, pixel size in width direction
    
    % get the height map
    input_table = readtable(file_name,'PreserveVariableNames',true);

    % pre-processing of the height map 
    table_x = size(input_table,2) - 1; % unit: num of pixels, -1: the last row is empty, table length
    matrix = table2array(input_table(:,2:table_x)); % convert the input_table to matrix and delete the last column.
    matrix = matrix/1000; % convert um to mm;
end

%% Function: to DownSample the Gel -> Match the Dimension of Raw Patterns of Cal dots
function hmap_ds = DownSampling(hmap,hmap_x,hmap_y,raw_x,raw_y)
    % creat an empty matrix to store the hmap after downsampling
    hmap_ds = zeros(raw_y, raw_x);

    % the size of the conversion matrix
    ratio_x = hmap_x/raw_x;
    ratio_y = hmap_y/raw_y;

    for i = 1:raw_y
        % define the start and end rows
        c_init = floor((i-1)*ratio_y + 1);
        c_term = ceil(i*ratio_y);
        
        if i == raw_y
            c_term = round(i*ratio_y);
        end

        % weight array for row
        [cur_size_c, cur_weight_c] = fitWeightArray(ratio_y,c_init,c_term,i);

        for j = 1:raw_x
            % define the start and end columns
            r_init = floor((j-1)*ratio_x + 1);
            r_term = ceil(j*ratio_x);
            
            if i == raw_y
                r_term = round(j*ratio_x);
            end

            % weight array for col
            [cur_size_r, cur_weight_r] = fitWeightArray(ratio_x,r_init,r_term,j);

            % current base matrix
            cur_grid = zeros(cur_size_c,cur_size_r);
            cur_grid(1:cur_size_c,1:cur_size_r) = hmap(c_init:c_term, r_init:r_term);

            % current weight matrix
            cur_weight = ones(cur_size_c,cur_size_r);
            cur_weight = cur_weight.*cur_weight_c';
            cur_weight = cur_weight.*cur_weight_r;

            % calculate avg
            cur_grid = cur_grid.*cur_weight;
            hmap_ds(i,j) = sum(cur_grid,'all')/sum(cur_weight,'all');
        end
    end
    
    % A Helper Function to get the Weight Array
    function [curr_size,current_weight] = fitWeightArray(ratio,init,term,idx)

        % obtain the current size of the weight array
        curr_size = term - init + 1;

        % weight of the start position  
        weight_str = ceil((idx-1)*ratio) - (idx-1)*ratio;
        if weight_str == 0
            current_weight(1,1) = 1;
        else
            current_weight(1,1) = weight_str;
        end
        
        % weight of the end position
        weight_end = idx*ratio - floor(idx*ratio);
        if weight_end == 0
            current_weight(end) = 1;
        else
            
            current_weight(end) = weight_end;
        end
    end
end

%% Function: to detrend the downsampling calibration gel height map
function hmap_dt = Detrending(hmap_ds)
    % get the dimension of hmap_ds
    y = size(hmap_ds,1);
    x = size(hmap_ds,2);

    % empty matrix to the t
    hmap_dt = zeros(size(hmap_ds));
    
    % calculate the length avg before detrending and the value need to be subtracted
    [x_avg_bef, x_dev] = fitModelOnAvgDt(hmap_ds,y,true);
    
    x_avg_aft = zeros(size(x_avg_bef));
    for i = 1:y 
        hmap_dt(i,:) = hmap_ds(i,:) - x_dev(i);
        x_avg_aft(1,i) = mean(hmap_dt(i,:));
    end
    
    % calculate the width avg before detrending and the value need to be subtracted
    [y_avg_bef, y_dev] = fitModelOnAvgDt(hmap_dt,x,false);

    y_avg_aft = zeros(size(y_avg_bef));
    for j = 1:x 
        hmap_dt(:,j) = hmap_dt(:,j) - y_dev(j);
        y_avg_aft(1,j) = mean(hmap_dt(:,j));
    end

    % helper function: Fit the Pattern with First Oreder Polynomial 
    function [avg_bef, deviation] = fitModelOnAvgDt(hmap_ds,oppo_size,is_length) % if is length, oppo_size should be width size
        
        % create an empty matrix to store the average value of each rol or col before detrending
        avg_bef = zeros(1,oppo_size);
        
        % calculate the average of row or col based on isRow value
        if is_length
            for len = 1:oppo_size
                avg_bef(len) = mean(hmap_ds(len,:));
            end
        else
            for col = 1:oppo_size
                avg_bef(col) = mean(hmap_ds(:,col));
            end
        end
        
        % fit two models for the calculated average
        dep_input = (1:oppo_size)';
        indep_input = avg_bef';
        fit_poly1 = fit(dep_input,indep_input,"poly1");
        deviation = polyval(coeffvalues(fit_poly1),1:oppo_size);
    end
end
 
%% Function: to perform cross correlation on different degree & different zooming scale
function [max_cor, deg_exp, sc_exp, offset, cor_values] = CalcCorrlation(gelH,gtH,mindeg,maxdeg,degdx,minsc,maxsc,scdx,chop)
    
    % set the range of the angle and the scaling coefficient
    scRange = minsc:scdx:maxsc;
    dgRange = mindeg:degdx:maxdeg;
    
    % empty matrix to store the cross correlation value
    cor_values = zeros(size(dgRange,2), size(scRange,2)); 
    xoffset_mat = zeros(size(dgRange,2), size(scRange,2));
    yoffset_mat = zeros(size(dgRange,2), size(scRange,2));
    
    % calculate cross correlation value at each degree and at each sacle coef
    for i = 1:size(dgRange,2)
        gelH_r = imrotate(gelH,dgRange(1,i),'bilinear');
        gelH_r_chop = gelH_r(chop+1:end-chop,chop+1:end-chop);
        for j = 1:size(scRange,2)
            % find the maximum correlation value
            gelH_s = imresize(gelH_r_chop,scRange(1,j)); 
            c = normxcorr2(gelH_s,gtH);
            % set up a boundary for maximum lag
            c_bound = c(size(gelH_s,1)+1:size(gtH,1),size(gelH_s,2)+1:size(gtH,2));
            cor_values(i,j) = max(c_bound,[],"all");
            % find the x-y coordinate of the max cross correlation value
            [peak_y_bound,peak_x_bound] = find(c_bound==max(c_bound(:)));
            peak_y = peak_y_bound + size(gelH_s,1);
            peak_x = peak_x_bound + size(gelH_s,2);
            % find the coordinate of the top left cornor of the gelsight
            xoffset = peak_x - size(gelH_s,2);
            yoffset = peak_y - size(gelH_s,1);
            xoffset_mat(i,j) = xoffset;
            yoffset_mat(i,j) = yoffset;
        end 
        disp(fprintf("degree complete %1.0f out of %1.0f", i, size(dgRange,2)))
    end
    
    max_cor = max(cor_values,[],"all");
    index = find(cor_values == max_cor); % index of the max correlation value
    offset = [xoffset_mat(index) yoffset_mat(index)];
    
    % find the scaling coef & rotation angles
    if index/size(dgRange,2) == floor(index/size(dgRange,2))
        sc_exp_idx = index/size(dgRange,2);
        deg_exp_idx = size(dgRange,2); 
    else
        sc_exp_idx = ceil(index/size(dgRange,2));
        deg_exp_idx = index - floor(index/size(dgRange,2))*size(dgRange,2); 
    end
    
    % output
    sc_exp = scRange(1,sc_exp_idx);
    deg_exp = dgRange(1,deg_exp_idx);
 end