%% Script Information

% Script Introduction:
%     This script create different square grating patterns for further usage. 

% Author:
%     Ronghao Zhang

% Dates:
%     Created: 2022-06-13
%     Last-Edited: 2022-08-24

% Methodology:
%    1. Create different square grating levels 
 
%% Initialization
% user defined parameters
window_size = 5; % (unit: mm) the diameter of the window area that contacts with the skin
touch_resolution = 0.025; % touch resolution (unit: mm) will be needed to create height maps
spat_freq = 1/touch_resolution; % unit: mm^(-1), spatial frequency per mm
pin_radius = touch_resolution/2; % unit:mm, the radiux of each 'pin' (smallest unit for sensation)

%% Defining the grating Strip Size

% define the basic parameters for the stripe
gap_width = [0.1, 0.5, 1.5, 3.0]; % gap_width (unit: mm) is defined as the width of the gap between gratings
grate_width = [0.1, 0.5, 1.0, 1.5];  % grate_width (unit:mm) is defined as the width of the grating
grate_height = [0.1, 0.2, 0.5, 1.0]; % grate_width (unit:mm) is defined as the height of the grating. 

% get the number of elements in all parameter array
num_gap_width = size(gap_width,2);
num_grate_width = size(grate_width,2);
num_grate_height = size(grate_height,2);
num_combination = num_gap_width*num_grate_height*num_grate_width;
grating_parameters = zeros(3,num_combination); 

for i = 1:num_gap_width
    % determine the mapping location of num_gap_width to the matrix
    idx_init = (i-1)*num_grate_width*num_grate_height + 1;
    idx_term = i*num_grate_width*num_grate_height;
    
    % create repeated values using elements in num_gap_width for mapping row 01 
    mapping_1 = repelem(gap_width(1,i), num_grate_width*num_grate_height); 
    
    % create repeated values using elements in num_grate_width for mapping row 02
    mapping_2 = repelem(grate_width,num_grate_height);
    
    % create repeated values using elements in num_grate_height for mapping row 03
    mapping_3 = repmat(grate_height, 1, num_grate_width);
    
    % assign the mapping vlaue to the matrix
    grating_parameters(1,idx_init:idx_term) = mapping_1(1,:);
    grating_parameters(2,idx_init:idx_term) = mapping_2(1,:);
    grating_parameters(3,idx_init:idx_term) = mapping_3(1,:);
end

%% Creating a Height Map for the grating Strip:

pin_size = window_size * spat_freq; % the number of each pins in one side 

height_map_cell = cell(1,num_combination); 

for i = 1:num_combination
    % get the grating strip information
    gap_width_pin = grating_parameters(1,i)*spat_freq; 
    grate_width_pin = grating_parameters(2,i)*spat_freq;
    grate_height = grating_parameters(3,i);
    
    % compute and store the strip information
    height_map = CalculateHeightMap(gap_width_pin, grate_width_pin, grate_height, pin_size);
    height_map_cell{1,i} = height_map;
end

%% Plot the Height Maps of Grating Patterns

figure; 
hold on
for i = 1:num_combination
    subplot(8,8,i)
    hold on
    hmap = height_map_cell{1,i};
    hmap_x = 0:touch_resolution:touch_resolution*size(hmap,2);
    hmap_y = 0:touch_resolution:touch_resolution*size(hmap,1);
    imagesc(hmap_x, hmap_y, hmap);
    axis([0 touch_resolution*size(hmap,1) 0 touch_resolution*size(hmap,2)]); 
    xlabel("width (mm)"); 
    ylabel("length (mm)");
    colormap hot;
    colorbar;
    caxis([0 1.1]);
    hold off
end
hold off

%% Save the height maps
save('C:\Users\somlab\Desktop\skin-biomechanics\data\others\square_gratings.mat', 'height_map_cell');

%% to create square gratings with differnt gap size, height, indent width
function hmap = CalculateHeightMap(dc, gw, gh, pz) % dc: discontinuity or gap width | gw: grate_width | gh: grate_height | wz: window_size | pz: pin_size
   
   hmap = ones(pz); % an empty height map
   
   num_set = floor(pz/(gw+dc)); % a set is one gap + one grate
   rmdr = (pz/(gw+dc) - num_set) * (gw+dc); % the reminder length of the window that is not covered by a complete set
   
   if gw >= dc % when the grate width >= the gap width -----------------------
       
       if rmdr >= gw % when the reminder larger than or equal to the grate width <<<
           
           %    |^^^^^^|    |^^^^^^|       |^^^^^^| 
           %    |      |    |      |   +   |      |    
           %----|      |----|      |       |      |--      
           
           grat_loc = zeros(1,pz);
           offset = ceil((rmdr-gw)/2);
           idx_init = offset + 1;
           for i = 1:(num_set+1)
               idx_term = idx_init + gw - 1; 
               grat_loc(1,idx_init:idx_term) = 1;
               idx_init = idx_term + 1 + dc; 
           end
           
       elseif rmdr > dc % when the reminder larger than the gap width <<<
           
           %    |^^^^^^|    |^^^^^^|       ^|
           %    |      |    |      |   +    |
           %----|      |----|      |        |----
           
           grat_loc = ones(1,pz);
           offset = ceil((rmdr-dc)/2); 
           idx_init = offset + 1;
           for i = 1:num_set
               idx_term = idx_init + dc - 1;
               grat_loc(1,idx_init:idx_term) = 0; 
               idx_init = idx_term + 1 + gw;
           end
           
       else % when the reminder smaller than both or equal to dc <<<
           
           %    |^^^^^^|    |^^^^^^|        |
           %    |      |    |      |   +    |
           %----|      |----|      |        |----
           
           grat_loc = ones(1,pz); 
           offset = floor((rmdr+gw)/2);
           idx_init = offset + 1; 
           for i = 1:num_set
              idx_term = idx_init + dc - 1;
              grat_loc(1,idx_init:idx_term) = 0; 
              idx_init = idx_term + 1 + gw;
           end
   
       end 
       
   else % when the gap width >= the grate width ------------------------------
       
       if rmdr > gw % when the reminder larger than the grate width <<< 
           
           %      |^^^^|      |^^^^|       |^^^^|
           %      |    |      |    |   +   |    |
           %------|    |------|    |       |    |---  
           
           grat_loc = zeros(1,pz);
           offset = ceil((rmdr-gw)/2);
           idx_init = offset + 1;
           for i = 1:(num_set+1)
               idx_term = idx_init + gw - 1; 
               grat_loc(1,idx_init:idx_term) = 1;
               idx_init = idx_term + 1 + dc; 
           end
           
       else % when the reminder smallwer than both <<< 
           
           %      |^^^^|      |^^^^|       
           %      |    |      |    |   +   
           %------|    |------|    |       ---  
           
           grat_loc = ones(1,pz);
           offset = ceil((rmdr+gw)/2);
           idx_init = offset + 1;
           for i = 1:num_set
               idx_term = idx_init + gw - 1; 
               grat_loc(1,idx_init:idx_term) = 0;
               idx_init = idx_term + 1 + gw; 
           end
       end
   end
   
   grat_hgt = grat_loc*gh;
   hmap = hmap.*grat_hgt; 
   
end




