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

    % empty matrix to store the detrended height map
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
        
        % fit a linear model
        dep_input = (1:oppo_size)';
        indep_input = avg_bef';
        fit_poly1 = fit(dep_input,indep_input,"poly1");
        deviation = polyval(coeffvalues(fit_poly1),1:oppo_size);
    end
end