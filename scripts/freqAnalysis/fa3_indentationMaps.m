%% Script Information

% Script Introduction:
%     This script aims to use the adjusted height map to generate
%     indentation plots and figures.

% Author:
%     Ronghao Zhang

% Dates:
%     Created: 2022-08-03
%     Last-Edited: 2022-08-05

% Methodology:
%     1. Load the preProcessing data that is previously generated and
%        stored using nci2_preProcessing script. 
%       *Remember to change the gel index. Options: M1/M3/M4
%     2. Plot the height maps 
%       *In this step, take the mean of the values > 95th percentile (do this for gel and ground truth)
%        and then add the difference to all gelsight map => this elevate the height of gel sight.
%     3. Plot the cross section profile for comparison. 
%     4. Save the modified gelsight height map (indentation map). 

%% Import the Preprocessed Image
load('C:\Users\somlab\Desktop\skin-biomechanics\data\freqAnalysis\squareGrating\m1_preProcessing.mat');

%% Create Color Themes
% N = 256; %size(get(gcf,'colormap'),1); % size of the current colormap
% vec = [      100;       75;       50;       25;        0];
% % Theme 01 -----
% hex1 = ['#0A0908';'#22333B';'#EAE0D5';'#F27059';'#994636'];
% raw1 = sscanf(hex1','#%2x%2x%2x',[3,size(hex1,1)]).' / 255;
% map1 = interp1(vec,raw1,linspace(100,0,N),'pchip');
% % Theme 02 -----
% hex2 = ['#2E4052';'#476C9B';'#EAE0D5';'#F27059';'#994636'];
% raw2 = sscanf(hex2','#%2x%2x%2x',[3,size(hex2,1)]).' / 255;
% map2 = interp1(vec,raw2,linspace(100,0,N),'pchip');
% % Theme 03 -----
% hex3 = ['#4A4A48';'#566246';'#EAE0D5';'#D66853';'#C75146'];
% raw3 = sscanf(hex3','#%2x%2x%2x',[3,size(hex3,1)]).' / 255;
% map3 = interp1(vec,raw3,linspace(100,0,N),'pchip');

% can call the colormaps using the command `colormap ...`

%% Plot the Height Map

m1_indent_match = cell(1,2);              %%%% CHANGE HERE M1/3/4

% title_str = ["M4 Ground Truth - Blue151"; 
%              "M4 Ground Truth - Chiffon";
%              "M4 Ground Truth - CrinkleySilk";
%              "M4 Gel Sight - Blue151";
%              "M4 Gel Sight - Chiffon";
%              "M4 Gel Sight - CrinkleySilk"];
         
% title_str = ["M1 Ground Truth - Denim"; 
%              "M1 Ground Truth - Microsuede";
%              "M1 Ground Truth - SandPaper";
%              "M1 Gel Sight - Denim";
%              "M1 Gel Sight - Microsuede";
%              "M1 Gel Sight - SandPaper"];

title_str = ["H0.25 W1.5 G1" "H0.25 W1.5 G3" "H0.25 W1.5 G5"];
% plot the ground truth and gel sight height map
figure;
hold on
for n = 1:3
    
    gt_hmap = m1_hmap_match{n,1};         %%%% CHANGE HERE M1/3/4
    gel_hmap = m1_hmap_match{n,2};        %%%% CHANGE HERE M1/3/4
    
    gt_max_avg = mean(gt_hmap(gt_hmap > prctile(gt_hmap,95,"all"))); 
    gel_max_avg = mean(gel_hmap(gel_hmap > prctile(gel_hmap,95,"all"))); 
    indent_hmap = gel_hmap + (gt_max_avg - gel_max_avg);
    
    m1_indent_match{n,1} = gt_hmap;       %%%% CHANGE HERE M1/3/4
    m1_indent_match{n,2} = indent_hmap;   %%%% CHANGE HERE M1/3/4
    
    subplot(2,3,n)
    hold on
    imagesc(gt_hmap);
    cmax = max(max(gt_hmap));
    cmin = min(min(gt_hmap));
    caxis([cmin cmax]);
    xlim([1 size(gt_hmap,2)]);
    ylim([1 size(gt_hmap,1)]);
    title(title_str(n));
    colormap hot
    colorbar
    hold off
    
    subplot(2,3,n+3)
    hold on
    imagesc(indent_hmap);
    caxis([cmin cmax]);
    xlim([1 size(gt_hmap,2)]);
    ylim([1 size(gt_hmap,1)]);
    title(title_str(n));
    colormap hot
    colorbar
    hold off
end
hold off

%% Plot the Cross Section of Gel and Grating  
% plot the ground truth and gel sight height map
figure;
hold on
for n = 1:3
    m1_gt_hmap = m1_indent_match{n,1};
    m1_indent_hmap = m1_indent_match{n,2};
    m3_gt_hmap = m3_indent_match{n,1};
    m3_indent_hmap = m3_indent_match{n,2};
    m4_gt_hmap = m1_indent_match{n,1};
    m4_indent_hmap = m1_indent_match{n,2};
    
    
    subplot(3,3,n)
    hold on 
    plot(1:size(m1_gt_hmap,1),m1_gt_hmap(:,1000),"red");
    plot(1:size(m1_gt_hmap,1),m1_indent_hmap(:,1000),"blue");
    ylim([0,3]);
    ylabel('Height (mm)');
    xlabel('Pixels');
    title('M1 Gel on the gratings with height 0.25 mm');
    hold off
    
    subplot(3,3,n+3)
    hold on 
    plot(1:size(m3_gt_hmap,1),m3_gt_hmap(:,1000),"red");
    plot(1:size(m3_gt_hmap,1),m3_indent_hmap(:,1000),"blue");
    ylim([0,3]);
    ylabel('Height (mm)');
    xlabel('Pixels');
    title('M3 Gel on the gratings with height 1 mm');
    hold off
    
    subplot(3,3,n+6)
    hold on 
    plot(1:size(m4_gt_hmap,1),m4_gt_hmap(:,1000),"red");
    plot(1:size(m4_gt_hmap,1),m4_indent_hmap(:,1000),"blue");
    ylim([0,3]);
    ylabel('Height (mm)');
    xlabel('Pixels'); 
    title('M4 Gel on the gratings with height 2.5 mm');
    hold off
    
end
hold off

%% Save the Indentation Map 
save('C:\Users\somlab\Desktop\nonCompliantTouchSim_rz_2022\data\squareGrating\combined_indent.mat', ...
     'm1_indent_match', 'm3_indent_match', 'm4_indent_match');
