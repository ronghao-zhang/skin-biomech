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
    
    % creat an empty matrix to storaq the hmap after downsampling
    hmap_ds = zeros(raw_y, raw_x);
    
    % calculate the current size in x and y dimension
    hmap_x = size(hmap,2);
    hmap_y = size(hmap,1);

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