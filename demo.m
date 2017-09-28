clear all; close all; clc;

addpath([pwd, '\dependencies\'])

%% Load data

% video stored in chunks in order to comply with GitHub file size limits
% load and concatenate video
loc   = [pwd, '\video_files\'];
N     = 1906;
video = zeros(300, 300, N, 'uint16');
len   = 191;
for ii = 1:10
   ids = (ii-1)*len + (1:len);
   ids = ids(ids <= N);
   load([loc, 'video_', num2str(ii)]);
   video(:,:,ids) = vid;
end

radius = 7;
meanIm = mean(video,3);      % mean image
corrIm = crossCorr(video);   % correlation image

%% Initialisation

retune_alpha  = 0;           % if 1, retune the value of alpha

% Parameters
alpha                       = 0.55;
init_opt.blur_radius        = 1.5;

if retune_alpha
    
    % Format figure that is used for tuning alpha
    init_opt.figureHandle = 1;
    figure;
    set(gcf, 'Position', [10, 10, 1910, 385]);      

    % Tune alpha for metric (correlation image)
    [alpha, phi_0] = tune_alpha(corrIm, radius, alpha, init_opt);

else
    phi_0          = initialise(corrIm, radius, alpha, init_opt);  
end

plotContoursOnSummaryImage(corrIm, phi_0<0, []);

%% Segmentation

retune_lambda = 0;

% Refresh workspace
close all; clear seg_opt;

% Algorithm parameters
seg_opt.lambda              = 10;
seg_opt.mergeCorr           = 0.95;
seg_opt.mergeDuring         = 1;

if retune_lambda  
	lambda         = tune_lambda(phi_0, video, radius,...
                                 seg_opt, corrIm, meanIm);
    seg_opt.lambda = lambda;
end

% Do segmentation, segmenting about 300 ROIs takes ~15 mins
tic;
[masks, cell_ts, nhbd_ts] = segment(phi_0, video, radius, seg_opt);
runtime = toc

% Get number of pixels in each mask
pix_num = zeros(size(masks,3),1);
for mask_num = 1:size(masks,3)
   pix_num(mask_num) = nnz(masks(:,:,mask_num));
end

%% Display initial results, before size thresholding

close all; clear opts

% masks on summary images
opts.subplot = 1;
opts.m       = 1;
opts.n       = 2;
opts.p       = 1;
plotContoursOnSummaryImage(corrIm, masks, opts);
title('Contours on correlation image');
opts.p       = 2;
% scale meanIm for better visibility
opts.scaled  = 1;
scaled_meanIm = meanIm - min(meanIm(:));
scaled_meanIm = 3.5*scaled_meanIm./max(scaled_meanIm(:));
plotContoursOnSummaryImage(scaled_meanIm, masks, opts);
title('Contours on mean image');

%% example of small and large ROI
T          = 1/8;
t          = T:T:(N*T);
line_width = 1.5;
small_ID   = 17;
large_ID   = 147;
figure;
opts.m = 2;
opts.n = 4;
opts.p = 1;
plotContoursOnSummaryImage(corrIm, masks(:,:,small_ID), opts);
axis([1, 40, 161, 200])
title('small ROI', 'FontSize', 18);
subplot(2, 4, 2:4)
plot(t, cell_ts(small_ID,:), 'LineWidth', line_width);
hold on
plot(t, nhbd_ts(small_ID,:), 'LineWidth', line_width);
xlabel('Time(s)', 'FontSize', 18);
set(gca, 'xtick', 0:50:200)
hdl = legend('Interior time series', 'Neighbourhood time series',...
       'Location', 'northoutside', 'Orientation', 'Horizontal');
set(hdl, 'FontSize', 18);
set(gca, 'ytick', [])
box off
legend boxoff
opts.p = 5;
plotContoursOnSummaryImage(corrIm, masks(:,:,large_ID), opts);
axis([161, 200, 21, 60])
title('large ROI', 'FontSize', 18);
subplot(2, 4, 6:8)
plot(t, cell_ts(large_ID,:), 'LineWidth', line_width);
hold on
plot(t, nhbd_ts(large_ID,:), 'LineWidth', line_width);
xlabel('Time(s)', 'FontSize', 18);
set(gca, 'xtick', 0:50:200)
set(gca, 'ytick', [])
hdl = legend('Interior time series', 'Neighbourhood time series',...
       'Location', 'northoutside', 'Orientation', 'Horizontal');
set(hdl, 'FontSize', 18);
box off
legend boxoff

%% Post processing (size thresholding)
figure;
% If a user is only looking for cell bodies, it is beneficial to threshold
% the size of the ROIs. 
min_size     = pi*radius^2*0.5;
smaller_ROIs = masks(:,:, pix_num < min_size);
larger_ROIs  = masks(:,:, pix_num >= min_size);

clear options
options.subplot = 1;
options.m       = 2;
options.n       = 2;
options.p       = 1;
h(1)            = plotContoursOnSummaryImage(corrIm, larger_ROIs, options);
title('larger ROIs on corrIm')
options.p       = 3;
h(3)            = plotContoursOnSummaryImage(corrIm, smaller_ROIs, options);
title('smaller ROIs on corrIm')
options.scaled  = 1;
options.p       = 2;
h(2)            = plotContoursOnSummaryImage(scaled_meanIm, larger_ROIs, options);
title('larger ROIs on meanIm')
options.p       = 4;
h(4)            = plotContoursOnSummaryImage(scaled_meanIm, smaller_ROIs, options);
title('smaller ROIs on meanIm')
linkaxes([h(1), h(2), h(3), h(4)])
