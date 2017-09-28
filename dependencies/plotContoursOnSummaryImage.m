function[hdl] = plotContoursOnSummaryImage(metric, cellMasks, opt)

% AUTHOR: Stephanie Reynolds (25/09/2017)
%
% OVERVIEW: This function plots the masks of the ROIs in cellMasks on the
% summary image ('metric').
%
%%%%%%%%%%%%%%%   INPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metric                  MxN summary image of video, usually the
%                         pixelwise cross-correlation (see crossCorr.m) or
%                         mean image
% cellMasks               MxNxK binary array, for the i-th sheet 
%                         (cellMasks(:,:,i)) a pixel's value is 1 if it is 
%                         in the i-th cell and it is 0 if that pixel is not
%                         in the i-th cell.
% opt                     A variable of type struct. In the following we
%                         describe the fields of this struct. 
% opt.subplot             If present, the figure is plotted as a subplot in
%                         location subplot(m, n, p) where opt.m = m, opt.n
%                         = n and opt.p = p.
% opt.plot_ids            If present, the ID of the ROI (the sheet of the
%                         cellMasks array to which it corresponds) is
%                         written as text next to the contour.
% opt.printPlot           If this field is present, the plot is printed in
%                         file location opt.figLoc with name opt.figName.
%
%%%%%%%%%%%%%%%   OUTPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hdl                    The figure handle.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cellMasks = logical(cellMasks);


if ~isfield(opt, 'subplot')
    hdl = figure;
elseif isfield(opt, 'figureHandle')
    gcf = opt.figureHandle;
else
    hdl = subplot(opt.m, opt.n, opt.p);
end

% Format figure
rng(1);
line_width   = 1.5;
cols         = mod(randperm(size(cellMasks,3)),5) + 1; 
auto_colour  = [161, 69, 10;...
                204, 97, 28;...
                230, 131, 66;...
                255, 165, 106;...
                255, 191, 150]./255;        


% Rescale and plot metric
if ~isfield(opt, 'scaled')
    metric = metric - min(metric(:));
    metric = metric/max(metric(:));
end
imshow(metric);
cmap = cmocean('tempo');
cmap = flip(cmap);
colormap(gca, cmap)
hold on

% Plot contours on metric
for ii = 1:size(cellMasks,3)
    
    colour     = auto_colour(cols(ii),:);
    [B,~]          = bwboundaries(cellMasks(:,:,ii));
    for jj = 1:length(B)
        boundary       = B{jj};
        plot(boundary(:,2), boundary(:,1),'Color', colour, 'LineWidth',... 
            line_width);
    end
    if isfield(opt, 'plot_ids')
       [y_loc, x_loc] = find(cellMasks(:,:,ii));
        text_y         = median(y_loc);
        text_x         = median(x_loc);
        text(text_x, text_y, num2str(ii), 'FontSize', 10, 'Color', colour);
        
    end
      
end
if ~isfield(opt, 'subplot')
    set(gca,'Position',[0 0 1 1]);
end
set(gca, 'xtick',[])
set(gca, 'ytick',[])
pbaspect([1,1,1]);

% Print
if isfield(opt, 'printPlot')
    if isfield(opt, 'figName')
        figName = opt.figName;
    else
        figName = 'contoursOnSummaryImage';
    end
    figName = [opt.figLoc, '\', figName]; 
    print(figName, '-dpng', '-r300');
    savefig([figName, '.fig']);
end


end