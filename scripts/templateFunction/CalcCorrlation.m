function [max_cor, deg_exp, sc_exp, offset, cor_values] = CalcCorrlation(gelH,gtH,mindeg,maxdeg,degdx,minsc,maxsc,scdx,chop)

    %------------------------------
    %-Function aim- 
    %   Calculate the maximum cross correlation value at different rotation
    %   angle and at different scaling coefficient. This function gives the
    %   best fit angle for rotation and the best fit scaling coefficent
    %   for zooming. 
    %-Function Input-
    %   gelH: the height map of the gel sight after being downsampled and detrended 
    %   gtH: the height map of the ground truth after being detrended
    %   mindeg: the minimum degree for rotating the gelsight height map
    %   maxdeg: the maximum degree for rotating the gelsight height map
    %           the range will be -maxdeg to maxdeg
    %   degdx: the step size of the degree
    %   minsc: the minimum scaling coefficint. Since gel "magnifies" the
    %          ground truth, the range of the scaling coef. should be 0-1.
    %          which means => gel sight * sc = ground truth. 
    %   maxsc: the maximum scaling coefficient.
    %   scdx: the step size of the scaling coeffienct. 
    %   chop: the number of pixels choped on the boundray of the
    %         rectangular height map. 
    %-Function Output-
    %   max_cor: the max cross correlation value
    %   deg_exp: the best fit correlation degree
    %   sc_exp: the best fit scaling coefficient
    %   offset: the coordinate of left upper corner of the gelsight when matching gelsight to ground truth
    %   cor_values: the cross correlation value matrix
    %------------------------------
    
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