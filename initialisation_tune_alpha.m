function[alpha, phi_0] = ...
    initialisation_tune_alpha(alpha, cellMetric, radius, options)
 
if alpha(1) == 0
    alpha(1) = 0.1;
end
if length(alpha) > 1
    if alpha(2) == 0
        alpha(2) = 0.1;
    end
end
 
temp_options              = options;
cellMetric                = cellMetric - min(cellMetric(:));
cellMetric                = cellMetric/max(cellMetric(:));

[initMasks] = initialise(cellMetric, radius, alpha, options);
while size(initMasks,3) >700 
    alpha       = 2*alpha;
    initMasks   = initialise(cellMetric, radius, alpha, options);
end
initMasks   = initMasks<0;
disp('Start tuning first alpha...');       
if isfield(options, 'figureHandle')
    figure(1);
    h1 = subplot(1,3,1);
    plotMetric = cellMetric;%imgaussfilt(cellMetric,2);
    imagesc(plotMetric)
else
    figure; 
    h1 = subplot(1,3,1);
    plotMetric = cellMetric;%imgaussfilt(cellMetric,2);
    imagesc(plotMetric)
    pause
end
pbaspect([1,1,1])
new_alpha     = alpha;
current_alpha = alpha;
new_alpha(1)  = 5/6*alpha(1);

h2 = subplot(1,3,2);
h3 = subplot(1,3,3);
linkaxes([h1,h2,h3])
while new_alpha(1) ~= current_alpha(1)
    subplot(1,3,2);
    imshowpair(plotMetric, any(initMasks,3))
    title(['alpha = ',num2str(current_alpha),...
           ', cell no = ', num2str(size(initMasks,3))])

    phi_0         = initialise(cellMetric, radius,...
                               new_alpha, temp_options);
    initMasks     = phi_0<0;
    subplot(1,3,3);
    imshowpair(plotMetric, any(initMasks,3))
    title(['new alpha = ',num2str(new_alpha(1)),...
           ', cell no = ', num2str(size(initMasks,3))])
    current_alpha = new_alpha;
    disp([num2str(size(initMasks,3)), ' masks initialised.']);
    new_alpha(1)  = input(['current alpha = ', num2str(current_alpha(1)),', what is new alpha? ']);
    

end
alpha(1) = new_alpha(1);    




phi_0 = initialise(cellMetric, radius, alpha, options);

end