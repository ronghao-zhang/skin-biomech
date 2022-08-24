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