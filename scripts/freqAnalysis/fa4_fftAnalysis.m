%% Script Information

% Script Introduction:
%     This script aims to perfrom fast fourier transformation, and to plot
%     the power density spectrum.

% Author:
%     Ronghao Zhang

% Dates:
%     Created: 2022-08-08
%     Last-Edited: 2022-08-24

% Check before Run
%   1. change the gel number
%   2. use Crtl + F to replace all m? with the correct gel index

%% Load the Data

m1_gel_num = 3; % number of pairs (1 pair = 1 gel + 1 non-gel)

% load the calibrated and paired height map
load('C:\Users\somlab\Desktop\skin-biomechanics\data\freqAnalysis\squareGrating\m1_preProcessing.mat'); 

%% Adjust Indentation Depth for GelSight (For Non-Compliant Object)

% empty array to store the hmap after ajust indentation.
m1_indent_match = cell(m1_gel_num,2); 

for n = 1:m1_gel_num
    
    % obtain ground truth and gelsight
    gt_hmap = m1_hmap_match{n,1};
    gel_hmap = m1_hmap_match{n,2};
    
    % calculate the average of values > 95 percentile for both gel and gt
    gt_max_avg = mean(gt_hmap(gt_hmap > prctile(gt_hmap,95,"all"))); 
    gel_max_avg = mean(gel_hmap(gel_hmap > prctile(gel_hmap,95,"all"))); 
    indent_hmap = gel_hmap + (gt_max_avg - gel_max_avg); % add the difference to current gel hmap
    
    % export 
    m1_indent_match{n,1} = gt_hmap;
    m1_indent_match{n,2} = indent_hmap;
end

% %% Perform Fast Fourier Transformation & Power Density Spectrum (Non-Compliant Object)
% 
% % Change the Directory to the Data Directory
% cd('C:\Users\somlab\Desktop\skin-biomechanics\raw\gelsight2022\squareGratings\m1\')
% 
% % Obtain Ground Truth information
% m1_grdtruth = dir(fullfile(pwd,'.\None*.csv'));
% 
% for n = 1%:m1_gel_num
%     
%     % obtain the processed patterns
%     gt_hmap = m1_indent_match{n,1};
%     indent_hmap = m1_indent_match{n,2};
%     
%     % obtain the resolution of the ground truth
%     f_name = m1_grdtruth(n).name; 
%     res = dlmread(f_name,',',[3,1,3,1])/1000; % unit: mm per pixel, pixel size in length direction
%     
%     [gt_freq_vec, gt_power] = CalcPowerPerFreq(gt_hmap,res);
%     [gs_freq_vec, gs_power] = CalcPowerPerFreq(indent_hmap,res);
% end
% 
% figure; 
% hold on
% plot(gt_freq_vec, 10*log10(gt_power), "blue");
% plot(gs_freq_vec, 10*log10(gs_power), "red");
% hold off
% 
% figure; 
% hold on
% plot(gt_freq_vec, gt_power, "blue");
% plot(gs_freq_vec, gs_power, "red");
% hold off
% 
% figure; 
% hold on
% plot(gt_freq_vec, gt_power/gs_power, "black");
% hold off
% 
% 
% %% Function: to Calculate Normalized Power Spectral Density for Square Grating Patterns
% 
% function [freq_vec, hmap_psd] = CalcPowerPerFreq(hmap,res)
%     % the sampling frequency
%     freq = 1/res;
%     
%     hmap_fft = mean(fftshift(fft(hmap)),2); % the average fft of the hmap 
%     hmap_psd = abs(hmap_fft).^2;
%     hmap_psd = hmap_psd./sum(hmap_psd(:));
%     freq_vec = 0:freq/size(hmap,1):freq;
%     freq_vec = freq_vec(1:end-1); 
% end
% 
% % correct code
% L = length(trace);
% trace_fft = fft(trace) * freq; % 1
% trace_fft_2sided = abs(trace_fft); % 2
% trace_fft_1sided = trace_fft_2sided(1:L/2+1); % 3
% trace_fft_1sided(2:end-1) = 2*trace_fft_1sided(2:end-1); % 4