function[alpha, phi_0] = tune_alpha(metric, radius, alpha, options)

% AUTHOR: Stephanie Reynolds (25/09/2017)
%
% OVERVIEW: This file enables a user to tune the value of alpha used in the
% initialisation algorithm. Tuning stops when the user-input value of alpha
% is the same as the existing value of alpha.
%
%%%%%%%%%%%%%%%   INPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha                 initial estimate of alpha, suggested value 0.5
% cellMetric            the MxN summary image (e.g. correlation image, 
%                       mean image,...) that is used for the initialisation
% radius                expected radius of a cell
% options               variable of type struct. see initialise.m for 
%                       further explanation
% options.figureHandle  if this field exists the plot will be done on the
%                       existing figure
%
%%%%%%%%%%%%%%%   OUTPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha       the value of the tuned parameter alpha
% phi_0       MxNxK array, where K is the number of initialised ROIs.
%             each sheet (phi_0(:,:,ii)) is -1 when a pixel is in the i-th ROI, 
%             and equal to 1 when a pixel is not in the i-th ROI

if alpha(1) == 0
    alpha(1) = 0.1;
end

% Rescale metric so it is more visible on plot
metric      = metric - min(metric(:));
metric      = metric/max(metric(:));

% Initialisation with first alpha
initMasks   = initialise(metric, radius, alpha, options);

% If we have lots of masks, increase alpha (=> more conservative)
while size(initMasks,3) >700 
    alpha       = 2*alpha;
    initMasks   = initialise(metric, radius, alpha, options);
end
initMasks   = initMasks<0; 

% Initialse plots
if isfield(options, 'figureHandle')
    h1 = subplot(1,3,1);
    imagesc(metric)
else
    figure; 
    h1 = subplot(1,3,1);
    imagesc(metric)
    disp('Resize figure window for best visibility.')
    pause
end
pbaspect([1,1,1])
opt.subplot   = 1;
opt.m         = 1;
opt.n         = 3;

% Start the tuning
disp('Start tuning first alpha...');      
new_alpha     = alpha;
current_alpha = alpha;
new_alpha(1)  = 5/6*alpha(1);
while new_alpha(1) ~= current_alpha(1)
    opt.p = 2;
    h2    = plotContoursOnSummaryImage(metric, initMasks, opt);
    title(['alpha = ',num2str(current_alpha),...
           ', cell no = ', num2str(size(initMasks,3))])

    phi_0      = initialise(metric, radius,...
                               new_alpha, options);
    initMasks  = phi_0<0;
    opt.p      = 3;
    h3         = plotContoursOnSummaryImage(metric, initMasks, opt);
    title(['new alpha = ',num2str(new_alpha(1)), ', cell no = ', num2str(size(initMasks,3))])
    current_alpha = new_alpha;
    linkaxes([h1,h2,h3])
    disp([num2str(size(initMasks,3)), ' masks initialised.']);
    new_alpha(1)  = input(['current alpha = ', num2str(current_alpha(1)),', what is new alpha? (if satisfied, input current alpha) ']);
    
    
end
alpha(1) = new_alpha(1);    
phi_0    = initialise(metric, radius, alpha, options);

end