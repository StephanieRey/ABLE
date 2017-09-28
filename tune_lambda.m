function[lambda] = tune_lambda(phi_0, video, radius,...
                               options, corrIm,...
                               meanIm)
                           
% AUTHOR: Stephanie Reynolds (25/09/2017)
%
% REFERENCE: Reynolds et al. (2016) ABLE: an activity-based level set 
% segmentation algorithm for two-photon calcium imaging data. eNeuro
%
% OVERVIEW: This function is used to tune the value of lambda on a small
% section of the dataset. 
%
%%%%%%%%%%%%%%%   INPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% algorithmPars           Structure of options for the segmentation
%                         algorithm.
% corrIm                  Correlation image, see crossCorr.m 
% meanIm                  Mean image.
% phi_0                   Initialisation: MxNxK array, where K is the
%                         number of ROIs found, see initialise.m . In each 
%                         sheet (masks(:,:,ii)) -1 indicates that a
%                         pixel is inside the i-th ROI, +1 indicates that
%                         the pixel is not inside. 
% video                   MxNxT video.
% radius                  Radius of a cell.
%
%%%%%%%%%%%%%%%   OUTPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lambda                  Value of tuned parameter lambda.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings for tuning lambda
options.corrIm        = corrIm;
options.meanIm        = meanIm;
options.plot_progress = 1;
options.maxIt         = 25;

% Select region on which to test segmentation
message = ['Draw a circle (of radius at least 3 times the radius ',...
           ' of one cell, around the ROI you would like to ',...
           ' test segmentation on.'];
disp(message);
hdl  = figure; 
mask = roipoly(corrIm);

% Get ROIs in the region
test_masks = [];
counter    = 1;
for ii = 1:size(phi_0,3)
    intersection = nnz(and(mask, phi_0(:,:,ii)<0));
    if intersection
        test_masks(:,:,counter) = phi_0(:,:,ii);
        counter                 = counter + 1;
    end
end

% Get ID of central ROI
opts.plot_ids = 1;
opts.figureHandle = hdl;
plotContoursOnSummaryImage(corrIm, test_masks<0, opts);
central_ID = input('What is the ID of the most central ROI? ');
options.cell_to_monitor = central_ID;

% Do segmentation
finito      = 0;
lambda(1)   = options.lambda;
counter     = 2;
msg = ['Press 1 if we should increase lambda (if there was little or no activity),\n',...
       'press 2 if we should decrease lambda (the evolution of the level set function - in the leftmost plot on the second row - was unstable) \n',...
       'press 0 to accept this lambda.\n '];
while ~finito
    segment(test_masks, video, radius, options);
    user_input = input(msg);
    if user_input == 1
        if length(lambda) == 1
            new_lambda = 1.5*lambda;
        elseif lambda(end) > lambda(end-1)
            new_lambda = 1.5*lambda(end);  
        else
            new_lambda = lambda(end) - (lambda(end) - lambda(end-1))./2;
        end
     elseif user_input == 2
        if length(lambda) == 1
            new_lambda = 0.5*lambda;
        elseif lambda(end) < lambda(end-1)
            new_lambda = 0.5*lambda(end);  
        else
            new_lambda = lambda(end) - (lambda(end) - lambda(end-1))./2;
        end
    else
        finito = 1;
    end
    lambda(counter)      = new_lambda;
    counter              = counter + 1;
    options.lambda       = new_lambda;
    close all;
end

lambda = options.lambda;