%% Script Information

% Script Introduction:
%     This script aims to use the calibration data from the height map of
%     gel and of the calibration dots. For each pair of gels and
%     calibration dots, the calibration need to be done once.

% Author:
%     Ronghao Zhang

% Dates:
%     Created: 2022-06-29
%     Last-Edited: 2022-08-24

% Methodology:
%    1. Update the downsampling method using a rectangular (not square) downsampling grid
%    2. Detrend the downsampled height map
%    3. Any Points below the 5th percentile is zero
%    4. Image Registration for matching the gel height map to the ground truth
%    5. Plot the offset of the gel map 
%    6. Fit and Validate the Linear Regression Model 

%% Import Calibration Data

% Change the Directory to the Data Directory
cd('C:\Users\somlab\Desktop\skin-biomechanics\raw\calibrationDots');

% Obtain the Information of Gel Calibration Data ".csv" Files
cal_gel = dir(fullfile(pwd,'M*.csv'));
gel_num = size(cal_gel,1);

% Obtain the Information of the Ground Truth Data (raw pattern of calibration dots)
cal_dot = dir(fullfile(pwd,'None*.csv'));

% Import the Height Map and Resolution of Raw Pattern of the Calibration Dot
[cal_dot_res_x, cal_dot_res_y, cal_dot_hmap] = ImportHeightMap(cal_dot.name, cal_dot.folder);
cal_dot_hmap_x = size(cal_dot_hmap,2);
cal_dot_hmap_y = size(cal_dot_hmap,1);

% Create Empty arrays to store the Gel Height Map & Import
Cal_gel_res_x = zeros(1,gel_num);
Cal_gel_res_y = zeros(1,gel_num);
Cal_gel_hmap = cell(1,gel_num);
Cal_gel_hmap_x = zeros(1,gel_num);
Cal_gel_hmap_y = zeros(1,gel_num);

for n = 1:gel_num
    % obtain file names and folder path
    f_name = cal_gel(n).name;
    f_path = cal_gel(n).folder;
    
    % Import the Gel Height Map
    [res_x, res_y, hmap] = ImportHeightMap(f_name, f_path);
    
    % Export the Resolution and Height Map
    Cal_gel_res_x(n) = res_x;
    Cal_gel_res_y(n) = res_y;
    Cal_gel_hmap{1,n} = flip(hmap,2);
    Cal_gel_hmap_x(n) = size(hmap,2);
    Cal_gel_hmap_y(n) = size(hmap,1);
end

%% Detrend the Ground Truth (Calibration Dots Height Map)
cal_dot_hmap_dt = Detrending(cal_dot_hmap); 

%% Pre-Processing of the Gel height Map (Downsampling and Detrending)

% Empty Cell Array to Store the downsampled height map
Cal_gel_hmap_dt = cell(1,gel_num);

% Subtract the 5th percentile value of the ground truth
cal_dot_5hgt = prctile(cal_dot_hmap_dt,5,"All");
cal_dot_hmap_td = cal_dot_hmap_dt - cal_dot_5hgt; 
cal_dot_hmap_td(cal_dot_hmap_td < 0) = 0; 

% DownSample the HeightMap
for n = 1: gel_num
    
    % get the input (hmap & dimension)
    hmap = Cal_gel_hmap{1,n};
    hmap_x = Cal_gel_hmap_x(1,n);
    hmap_y = Cal_gel_hmap_y(1,n);
    
    % calculate the downsampled height map
    hmap_ds = DownSampling(hmap,cal_dot_hmap_x,cal_dot_hmap_y);
    
    % calculate the detrended height map
    hmap_dt = Detrending(hmap_ds); 
    
    % calculate the pre-indent
    hmap_dt_5hgt = prctile(hmap_dt,5,"All");
    hmap_dt_td = hmap_dt - hmap_dt_5hgt;
    hmap_dt_td(hmap_dt_td < 0) = 0;
    
    % export the detrended height map of the gel
    Cal_gel_hmap_dt{1,n} = hmap_dt_td; 
end

clear hmap hmap_x hmap_y hmap_ds hmap_dt hmap_dt_5hgt hmap_dt_td

%% Image Registration of the Downsampled And Detrended Height Maps

points = 40; % points that need to be manually select 

% Transpose the Dot Hmap Matrix
fix_hmap = cal_dot_hmap_td';
fix_points = SelectDotsOnFig(points, fix_hmap, 'Select 40 points in the Ground Truth Image. REMEMBER THE ORDER!!!');

img_size_x = size(fix_hmap,1); 
img_size_y = size(fix_hmap,2);

% Empty Cell to Store the Transformed Gel Matrix
Gel_hmap_trans = cell(1,gel_num);

% Empty Cell to Store the Coordinates of the Moving Points and Fixed Points
Points_coord = cell(1,gel_num); 

for n = 1%:gel_num
    
    % transpose the gel matrix
    mov_hmap = Cal_gel_hmap_dt{1,n}'; 
    mov_points = SelectDotsOnFig(points, mov_hmap, 'Determine Center on M1/M3/M4/M5 GelSight Height Map. NEED TO BE IN THE SAME ORDER AS GROUND TRUTH !!!');
    
    
    % transformation matrix
    trans_matrix = fitgeotrans(mov_points, fix_points, "nonreflectivesimilarity");
    
    % recovered height map
    Roriginal = imref2d([img_size_x, img_size_y]);
    rec_hmap = imwarp(mov_hmap,trans_matrix,'OutputView',Roriginal);
    
    % store the recovered height map
    Gel_hmap_trans{1,n} = rec_hmap;
    
    % store the Coordinates of points 
    Points_coord{1,n} = mov_points;
end

clear Roriginal rec_hmap trans_matrix

%% Gel Offset Calculation 

% Empty Cell Array to Store the Height of the Dots
Point_hgt = cell(2,gel_num);

% Empty cell Array to Store the Slope and Intercept
Slope_intcpt = cell(1,gel_num);
Mdl = cell(1,gel_num);

for n = 1:gel_num
    
    % moving points and fixed points for the height map
    mov_points = Points_coord{1,n};  
    
    % round the value for indexing
    fix_idx = round(fix_points);
    mov_idx = round(mov_points);
    
    % empty array to store the height of the point in hmap
    gel_point_hgt = zeros(points,1);
    dot_point_hgt = zeros(points,1);
    
    % get the height of each point
    for p = 1:points 
        dot_point_hgt(p,1) = fix_hmap(fix_idx(p,2), fix_idx(p,1));
        gel_point_hgt(p,1) = mov_hmap(mov_idx(p,2), mov_idx(p,1));
    end
    
    % fit a linear regression model
    indept = dot_point_hgt;
    respon = gel_point_hgt;
    mdl = fitlm(indept, respon);
    slope_intcpt = table2array(mdl.Coefficients(1:2,1));
    
    % store the slope and intercept
    Slope_intcpt{1,n} = slope_intcpt; 
    Mdl{1,n} = mdl;
    
    % store the height of the dots
    Point_hgt{1,n} = dot_point_hgt;
    Point_hgt{2,n} = gel_point_hgt; 
    
end

clear mov_points fix_idx mov_idx gel_point_hgt dot_point_hgt indept respon slope_intcpt mdl

%% Save Data in the Data Folder
save('C:\Users\somlab\Desktop\skin-biomechanics\data\freqAnalysis\calibration\cal_linear_models.mat', 'Slope_intcpt');

%% Calibration Dot Height Plot with Linear Regression Model
figure;
hold on
for n = 1:gel_num
    subplot(1,gel_num,n)
    hold on
    scatter(Point_hgt{1,n}, Point_hgt{2,n});
    plot(0:1:1.5, 0:1:1.5);
    plot(Mdl{1,n});
    axis([0 1.5 0 1]);
    xlabel("height in the raw pattern (mm)");
    ylabel("height in the gel (mm)");
    hold off
end
hold off

%% Function: to import the resolution and Height Map
function [res_x, res_y, matrix] = ImportHeightMap(file_name,file_path)

    %------------------------------
    %-Function aim- 
    %   Import Height Map Generated Using a Confocal Microscopy Scanning
    %   Over Some Texture. 
    %-Function Input-
    %   file_name: the name of the .csv file generated after scanning
    %   file_path: the path of the .csv file generated after scanning
    %-Function Output-
    %   res_x: the resolution of the height map in the x-axis direction (length)
    %          unit: um/pixel
    %   res_y: the resolution of the height map in the y-axis direction (width)
    %          unit: um/pixel
    %   matrix: the matrix that stores the height of the texture. The value
    %           of each element is the height at that specific location. 
    %------------------------------

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
function hmap_ds = DownSampling(hmap,raw_x,raw_y)
    
    %------------------------------
    %-Function aim- 
    %   Downsample the height map (name: A) of one texture (e.g. Gelsight) to match the
    %   dimension of another height map (name: B) (e.g. Ground Truth). The resulting
    %   height map A' has exactly the same dimension as B. 
    %-Function Input-
    %   hmap: the height map need to be downsampled (A).
    %   raw_x: the x dimension of the height map being referenced (B).
    %          can be computed by size(B,2).
    %   raw_y: the y dimension of the height map being referenced (A).
    %          can be computed by size(B,1).
    %-Function Output-
    %   hmap_ds: the height map after downsampling (A'). 
    %------------------------------
    
    % calculate the current size in x and y dimension
    hmap_x = size(hmap,2);
    hmap_y = size(hmap,1);
    
    % creat an empty matrix to store the hmap after downsampling
    hmap_ds = zeros(raw_y, raw_x);

    % the size of the conversion matrix
    ratio_x = hmap_x/raw_x;
    ratio_y = hmap_y/raw_y;

    for i = 1:raw_y
        % define the start and end rows
        c_init = floor((i-1)*ratio_y + 1);
        c_term = ceil(i*ratio_y);

        % weight array for row
        [cur_size_c, cur_weight_c] = fitWeightArray(ratio_y,c_init,c_term,i);

        for j = 1:raw_x
            % define the start and end columns
            r_init = floor((j-1)*ratio_x + 1);
            r_term = ceil(j*ratio_x);

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

    
    %------------------------------
    %-Function aim- 
    %   Detrend a height map to subtract a best-fit linear model from your
    %   data. This step makes one to more focus on the trend of the pattern
    %   itself instead of artifacts being introduced. 
    %-Function Input-
    %   hmap_ds: the height map need to be detrended. 
    %            usually, the input height map is already downsampled. 
    %-Function Output-
    %   hmap_dt: the height map after detrending. 
    %-Note-
    % A Pre-Indentation need to be added to all elements 
    % or the mean of the hmap_dt will be around 0.  
    %------------------------------
    
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

%% Function: To Select Points from A Height Map
function points_coord = SelectDotsOnFig(points, hmap, titleStr)

    %------------------------------
    %-Function Input-
    %   points: number of points need to be select
    %   hmap: the height map
    %   titleStr: the title of the figure
    %-Function Output-
    %   points_coord: the coordinates of the points
    %------------------------------

    % empty matrix to store the points selected
    points_coord = zeros(points,2);
    
    % get the size of the height map
    hmap_x = size(hmap,1); 
    hmap_y = size(hmap,2);
    
    % plot the height map for point selection
    figure('Position',get(0,'ScreenSize')); 
    hold on
    imagesc(hmap); 
    axis([0 hmap_y 0 hmap_x]);
    set(gca, 'PlotBoxAspectRatio', [1,1,1]);
    title(sprintf('%s (Choose %d Points)',titleStr, points));
    hold off
    
    % point selection and visualization
    for p = 1:points
        [x_temp, y_temp] = ginput(1);
        hold on
        scatter(x_temp, y_temp, 'filled'); 
        points_coord(p,:) = [x_temp, y_temp];
    end
end 
