function[masks] = initialise(metric, radius, alpha, options)

% AUTHOR: Stephanie Reynolds (25/09/2017)
%
% REFERENCE: Reynolds et al. (2016) ABLE: an activity-based level set 
% segmentation algorithm for two-photon calcium imaging data. eNeuro
%
% OVERVIEW: This function is used as the initialisation for the segmentation
% algorithm. Peaks in the 2D summary image(s) are identified as candidate
% ROIs. Peaks are found with an built-in MATLAB image processing function 
% 'imextendedmax'. Peaks with relative height (with respect to their
% neighbours) that are less than alpha x sigma are suppressed. Here, sigma is
% the standard deviation of the summary image and alpha is a tuning 
% parameter. 
%
%%%%%%%%%%%%%%%   INPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metric                     MxN summary image of video, usually the
%                            pixelwise cross-correlation (see crossCorr.m)
% radius                     radius of a cell
% alpha                      tuning parameter, peaks below alpha*sigma will be
%                            suppressed
% options                    A variable of type struct. In the following we
%                            describe the fields of this struct. 
% options.blur_radius        [Default: 1] If present, this is the radius of
%                            blurring applied to the summary images. 
%                            blurred with radius options.blur_radius.
% options.secondary_metric   M x N array corresponding to a summary image,
%                            e.g. the mean image. If this field is present, 
%                            initialisation is  performed on both the first 
%                            argument 'metric' and a second summary image. The value
%                            options.second_alpha (the value of alpha for the
%                            second summary image) must also be present. 
%
%%%%%%%%%%%%%%%   OUTPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% masks                      MxNxK array, where K is the number of ROIs found.
%                            In each sheet (masks(:,:,ii)): -1 indicates that a
%                            pixel is inside the i-th ROI, +1 indicates that the
%                            pixel is outside the i-th ROI.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim     = size(metric);
maxSize = round(pi * radius^2 * 1.5);

% Blur the input image (so that only peaks wider than 1 pixel are detected)
if isfield(options, 'blur_radius')
    blur_radius =  options.blur_radius;
else 
    blur_radius = 1;
end
metric          = imgaussfilt(metric, blur_radius);

% Find peaks of metric
h                     = sqrt(nanvar(metric(:))); 
metric(isnan(metric)) = 0;
BW                    = imextendedmax(metric, alpha*h);

% Find peaks of second metric, if required
if isfield(options, 'secondary_metric')
    secondaryMetric               = options.secondary_metric;
    secondaryMetric               = imgaussfilt(secondaryMetric, blur_radius);
    h                             = sqrt(nanvar(secondaryMetric(:))); 
    secondaryMetric(isnan(secondaryMetric)) = 0;
    alpha_2                       = options.second_alpha;
    BW                            = or(BW, imextendedmax(secondaryMetric, alpha_2*h));
end

% Each 4-connected component is an ROI
CC         = bwconncomp(BW, 4);
obj_num    = CC.NumObjects;
pixels     = CC.PixelIdxList;
masks      = zeros(dim(1), dim(2), obj_num);
for ii  = 1:CC.NumObjects
    mask             = zeros(dim);
    mask(pixels{ii}) = 1;
    masks(:,:,ii)    = mask;
end

% Remove any that are too large 
nnzs                          = squeeze(sum(sum(masks,1),2));
masks(:,:,nnzs > maxSize)     = [];

% Final ROIs
masks                         = -1*masks + ~masks;

end









