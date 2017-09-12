function[masks] = initialise(cellMetric, radius, alpha, options)

%%%%% Initialise ROIs for algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The initialisation procedure finds the peaks in the 2D summary image(s)
% to identity ROIs in the video. This is done with an built-in fast MATLAB
% image processing function 'imextendedmax'. We suppress peaks with
% relative height above their neighbours less than alphaxsigma, where sigma
% is the standard deviation of the summary image and alpha is a tuning 
% parameter. We subsequently merge any sufficiently close and correlated 
% ROIs, remove large ROIs (which tend to correspond to intensity 
% inhomogeneities) and expand small ROIs. 
% INPUTS
%
% cellMetric             :MxN summary image of video, usually pick avergae 
%                         pixel cross-correlation - this metric shows up
%                         cells well.
% radius                 :expected radius of a cell
% alpha                  :tuning parameter, peaks below alpha*sigma will be
%                         suppressed.

% OUPUTS
%
% masks                  :MxNxC array, where C is the number of ROIs found.
%                         In each sheet (masks(:,:,ii)) -1 indicates that a
%                         pixel is inside the ROI, +1 indicates that the
%                         pixel is outside.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(options, 'expand_small')
    expand_small = options.expand_small;
else
    expand_small = 0;
end

if isfield(options, 'blur')
     cellMetric  = imgaussfilt(cellMetric, options.blur_radius);
end
dim             = size(cellMetric);
if isfield(options, 'minSize')
    minSize = options.minSize;
else
    minSize         = round(pi * radius^2 * 0.25);
end
if isfield(options, 'maxSize')
    maxSize = options.maxSize;
else
    maxSize         = round(pi * radius^2 * 1.5);
end


%%%% Find peaks of cell metric
h       = sqrt(nanvar(cellMetric(:))); 
cellMetric(isnan(cellMetric)) = 0;
BW      = imextendedmax(cellMetric, alpha*h);

if isfield(options, 'second_metric')
    bkgdMetric                    = options.bkgd_metric;
    h          = sqrt(nanvar(bkgdMetric(:))); 
    bkgdMetric(isnan(bkgdMetric)) = 0;
    alpha_2    = options.second_alpha;
    BW         = or(BW, imextendedmax(bkgdMetric, alpha_2*h));
end

%%%% Each 4-connected component is a contender region
CC         = bwconncomp(BW, 4);
obj_num    = CC.NumObjects;
pixels     = CC.PixelIdxList;
masks      = zeros(dim(1), dim(2), obj_num);
for ii  = 1:CC.NumObjects
    mask             = zeros(dim);
    mask(pixels{ii}) = 1;
    masks(:,:,ii)    = mask;
end

%%% Remove any that are too large (these are likely to correspond to
%%% intensity inhomogeneities).
nnzs                          = squeeze(sum(sum(masks,1),2));
masks(:,:,nnzs > maxSize)     = [];
nnzs(nnzs > maxSize)          = [];
obj_num                       = size(masks,3);

if expand_small
    %%%% Expand any that are too small
    d        = options.expand_radius;
    se       = strel('disk', d);
    for ii = 1:obj_num
        if nnzs(ii) < minSize
            masks(:,:,ii) = imdilate(masks(:,:,ii), se);
        end
    end
end

masks = -1*masks + ~masks;

end









