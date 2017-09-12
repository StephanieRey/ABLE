function[cellMasks, timeSeries] = segment(phi_0, video, radius, algorithmPars)
 
% Code written by Stephanie Reynolds 12/09/2016.

%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% phi_0         : M x N x K initialsation, where M x N is the size of one 
%                 frame of the video.  Each sheet phi_0(:,:,ii) is the  
%                 initialisation of one ROI. phi_0(x,y,ii) < 0  if pixel 
%                 (x,y) is in the interior of ROI ii, and phi_0(x,y,ii) > 0  
%                 if pixel (x,y) is in the exterior.
% video         : M x N x T video to be segmented.
% radius        : Radius of a cell (integer, units are pixels).
% algorithmPars : Structure defining optional algorithm parameters. 
                      
%%%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cellMasks     : M x N x K' array. Each sheet cellMasks(:,:,ii) is the  
%                 binary mask of one detected ROI. 
% timeSeries    : average time course of each detected ROI.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% ALGORITHM PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if isfield(algorithmPars, 'delta')
    delta  = algorithmPars.delta;
else
    delta  = 1;
end
if isfield(algorithmPars, 'bandWidth')
    bandWidth  = algorithmPars.bandWidth;
else
    bandWidth  = 2 * radius;
end
if isfield(algorithmPars, 'timestep')
    timestep  = algorithmPars.timestep;
else
    timestep  = 10;
end
if isfield(algorithmPars, 'lambda')
    lambda = algorithmPars.lambda / timestep;
else
    lambda = 100 / timestep;
end
if isfield(algorithmPars, 'minimumSize')
    minimumSize = algorithmPars.minimumSize;
else
    minimumSize  = 3; 
end 
if isfield(algorithmPars, 'maximumSize')
    maximumSize = algorithmPars.maximumSize;
else
    maximumSize  = round(pi * radius^2 * 4); 
end 
if isfield(algorithmPars, 'noChangeNum')
    noChangeNum      = algorithmPars.noChangeNum;
else
    noChangeNum      = 40;
end
if isfield(algorithmPars, 'overlapSaturation')
    overlapSaturation = algorithmPars.overlapSaturation;
else
    overlapSaturation = 0;
end
if isfield(algorithmPars, 'mergeCorr')
    mergeCorr = algorithmPars.mergeCorr;
elseif isfield(algorithmPars, 'snr')
    mergeCorr = 0.8/(1+(10^(-snr/10)));
else
    mergeCorr = 0.8;
end
if isfield(algorithmPars, 'notifications_on')
    notifications_on  = algorithmPars.notifications_on;
else
    notifications_on  = 0;
end
if isfield(algorithmPars, 'maxIt')
    maxIt  = algorithmPars.maxIt;
else
    maxIt  = 150;
end 

%%%% Fixed algorithm parameters    
mu            = 0.2/timestep;
c0            = 3;
epsilon       = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
%%%%%%%% Initialise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

it             = 1;
t_len          = size(video,3);
cell_num       = size(phi_0,3);
active         = true(cell_num,1);
activity       = inf*ones(cell_num,noChangeNum);
cellIndex      = 1:cell_num;
phi            = phi_0;
video_dim      = [size(video,1), size(video,2)];
timeSeries     = zeros(cell_num, t_len, 'single');
region         = zeros([video_dim, cell_num]);
se_narrowband  = strel('square', round(bandWidth));
se_overlap     = strel('square', ceil(radius/2));
video_reshaped = reshape(video,[video_dim(1)*video_dim(2)*t_len,1,1]);
mods           = (0:1:(t_len - 1))*video_dim(1)*video_dim(2);
          
%%%%%%%%%   Initialise time series, phi and narrowband of each cell   %%%%%
ii = 1;
while ii <= size(phi,3)
     
   currentPhi          = phi(:,:,ii);
   currentMask         = currentPhi<0;
   timeSeries(ii,:)    = extractTimeSeries(currentMask, mods, video_reshaped, t_len);
   phi(:,:,ii)         = maskToSignedDistanceFunction(currentPhi, c0);
    
   if ~any(currentMask)
        phi(:,:,ii)            = [];
        timeSeries(ii,:)       = [];
        region(:,:,ii)         = [];
        active(ii)             = [];
        activity(ii, :)        = [];
       
        cellIndex(ii)          = [];
        if notifications_on
            disp(['Removed ROI ',num2str(ii), ' too similar too background.']); 
        end
   else
       ii                   = ii + 1;
   end
    
end
 
%%% If we're plotting progress, initialise the plots
if isfield(algorithmPars, 'plot_progress')
    figure;
    if isfield(algorithmPars, 'plot_multiple')
        plot_multiple  = algorithmPars.plot_multiple;
    else
        plot_multiple  = 1;
    end
    corrIm = algorithmPars.corrIm;
    meanIm = algorithmPars.meanIm;
    if isfield(algorithmPars, 'cell_to_monitor')
        cell_to_monitor  = algorithmPars.cell_to_monitor;
    else
        cell_to_monitor  = 1;
    end
    ncols = 5;
    nrows = 4;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% UPDATE: While iterations less than maximum and cells active %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
while and(it < maxIt, any(active))
     
    phi   = NeumannBoundCond(phi, video_dim); % Make phi satisfy neumann boundary conditions
    ii    = 1;
    
    % Are any cells touching and v. correlated? If yes, merge them
    if and(isfield(algorithmPars, 'mergeDuring'), length(cellIndex) > 1)
        correlation       = corrcoef(timeSeries');
        correlation       = tril(correlation, -1); % get rid of duplicates (it is symmetric)
        [p, q]            = find(correlation > mergeCorr); 
        while and(~isempty(p), ~isempty(q))
            % if they are overlapping (use square shape dilation
            % ebcause its faster)
            overlap = and( imdilate(phi(:,:,p(1))<0,se_overlap), ...
                           imdilate(phi(:,:,q(1))<0,se_overlap));
            if any(overlap(:))
                ordering      = sort([p(1),q(1)]); % we merge to the first one
                merge_remove = ordering(2); 
                merge_keep   = ordering(1);
              
                [active, activity, phi, region,...
                 timeSeries, cellIndex] =...
                 merge(merge_keep, merge_remove, phi, c0, timeSeries,...
                       video_reshaped, t_len,...
                       mods, se_narrowband,...
                       active, activity,...
                       region, cellIndex); 

                p(1)               = [];
                q(1)               = [];
                p(p==merge_remove) = merge_keep;
                q(q==merge_remove) = merge_keep;
                p(p>merge_remove)  = p(p>merge_remove) - 1;
                q(q>merge_remove)  = q(q>merge_remove) - 1;
                unique_pq          = unique([p, q], 'rows');
                if ~isempty(unique_pq)
                    p                  = unique_pq(:,1);
                    q                  = unique_pq(:,2);
                end

            else
                p(1) = [];
                q(1) = [];
            end
        end    
    end
        
    while ii <= size(phi,3)
         
       if active(ii)  
           
           phi_removed    = 0;
           phiUpdate      = phi(:,:,ii);
           nhbd           = imdilate(phiUpdate<0, se_narrowband);
           region(:,:,ii) = nhbd;
           otherCells     = zeros(video_dim);
           ind            = find(nhbd);
           [ind_x, ind_y] = ind2sub(video_dim,ind) ;
           for jj = 1:length(ind_x)
               otherCells(ind_x(jj), ind_y(jj)) = any(phi(ind_x(jj),ind_y(jj),1:end ~= ii)<0,3);
           end
           nhbd             = and(nhbd, ~or(otherCells, phiUpdate<0));
           n                = extractTimeSeries(nhbd, mods, video_reshaped, t_len);
           diracPhi         = Dirac(phiUpdate, epsilon); 
           evaluate_idx     = diracPhi > 0;
           onlyThisMask     = and(~otherCells, phiUpdate<0);
           timeSeries(ii,:) = extractTimeSeries(onlyThisMask, mods,...
                                                video_reshaped, t_len);
           
            
           if nnz(evaluate_idx)
                [similarityVelocity] = ...
                cellSimilarityVelocity(video_dim, phi, n, ii, evaluate_idx,...
                                       mods, video_reshaped,  t_len, ...
                                       timeSeries, overlapSaturation,...
                                       otherCells);
           else
                similarityVelocity = zeros(size(diracPhi));
           end
           
           %%% Calculate regularizer
           regularisationVelocity = zeros(video_dim);
           [x_loc, y_loc]         = find(region(:,:,ii));
           x_min                  = max(min(x_loc) - 2, 1);
           x_max                  = min(max(x_loc) + 2, video_dim(2));
           y_min                  = max(min(y_loc) - 2, 1);
           y_max                  = min(max(y_loc) + 2, video_dim(2));
           if and(x_max > x_min, y_max > y_min)
              regularisationVelocity(x_min:x_max, y_min:y_max) =...
                distanceRegularisation(phiUpdate(x_min:x_max, y_min:y_max));
           end
          
           %%% Update phi
           phiUpdate              = phiUpdate + ...
                                    timestep * region(:,:,ii).*...
                                    (mu * regularisationVelocity - ...
                                    lambda*diracPhi.*...
                                    similarityVelocity);
                                
           degenerate   = ~any(phiUpdate(:)<0);
           onlyThisMask = and(~otherCells, phiUpdate<0);
           mask_size    = nnz(phi(:,:,ii)<0);                            
 
           
           if and(and(~degenerate, any(onlyThisMask(:))), mask_size < maximumSize)
 
               %%% Update stored values
               changed            = nnz(or(and(phi(:,:,ii)<0, phiUpdate>0),...
                                           and(phi(:,:,ii)>0, phiUpdate<0)));
               activity(ii,:)     = [activity(ii,2:end), changed];
               phi(:,:,ii)        = phiUpdate; 
 
               %%% If it has converged or is too big, then STOP (but don't remove).                                          
               if and(it > noChangeNum, activity(ii,:) < delta)
                    if mask_size < minimumSize
                        cellIndex(ii)         = [];
                        phi_removed           = 1;
                        phi(:,:,ii)           = [];
                        timeSeries(ii,:)      = [];
                        active(ii)            = [];
                        activity(ii, :)       = [];
                        region(:,:,ii)        = [];
                        if notifications_on
                            disp(['Removed ROI ',num2str(ii),' final size too small']);
                        end
                    else
                        active(ii)       = 0;
                        timeSeries(ii,:) = extractTimeSeries(onlyThisMask, mods,video_reshaped, t_len);
                        
                        if notifications_on
                            disp(['ROI ', num2str(ii), ' converged.']);
                        end
                    end
               end
           else %%% Is update degenerate? Does this mask entirely overlap with another or is it too too small? If so remove.
               %%% Display why ROI was removed
               if notifications_on
                   if ~any(onlyThisMask(:))
                      disp(['Removed ROI ',num2str(ii), ' completely overlapping.']);
                   elseif mask_size > maximumSize
                       disp(['Removed ROI ',num2str(ii), ' outsized.']);
                   else
                       disp(['Removed ROI ',num2str(ii), ' degenerate.']);
                   end
               end
               cellIndex(ii)        = [];
               phi_removed          = 1;
               phi(:,:,ii)          = [];
               timeSeries(ii,:)     = [];
               region(:,:,ii)       = [];
               active(ii)           = [];
               activity(ii, :)      = [];
           end
           if isfield(algorithmPars, 'plot_progress')
               if and(mod(it, plot_multiple) == (plot_multiple-1), cellIndex(ii) == cell_to_monitor)
                      updatePlot(ii, video, phi, radius, nrows, ncols, corrIm,...
                         meanIm, onlyThisMask, nhbd, region, diracPhi,...
                         similarityVelocity, regularisationVelocity,...
                         mu, mods, video_reshaped, t_len, timeSeries,...
                         it, activity, delta, noChangeNum)   
                      drawnow;
               end
           end 
                     
           if ~phi_removed
               ii = ii + 1;
           end
           
       else
           %%%%% Update time series if any cell has moved into it
           otherCells     = zeros(video_dim);
           [ind_x, ind_y] = find(phi(:,:,ii)<0);
           for jj = 1:length(ind_x)
               otherCells(ind_x(jj), ind_y(jj)) = any(phi(ind_x(jj),ind_y(jj),1:end ~= ii)<0,3);
           end
           onlyThisMask     = and(~otherCells, phi(:,:,ii)<0);
           timeSeries(ii,:) = extractTimeSeries(onlyThisMask, mods,video_reshaped, t_len);
           ii               = ii + 1;
       end

    end
    %toc;
    if notifications_on
        disp(['Iteration ', num2str(it)]);
    end
    it = it + 1; 
     
end

it = it - 1;
 
% Are any cells touching and v. correlated? If yes, merge them
if and(isfield(algorithmPars, 'mergeAtEnd'), length(cellIndex) > 1)
    correlation       = corrcoef(timeSeries');
    correlation       = tril(correlation, -1); % get rid of duplicates (it is symmetric)
    [p, q]            = find(correlation > mergeCorr); 
    while and(~isempty(p), ~isempty(q))
        % if they are overlapping (use square shape dilation
        % ebcause its faster)
        overlap = and( imdilate(phi(:,:,p(1))<0,se_overlap), ...
                       imdilate(phi(:,:,q(1))<0,se_overlap));
        if any(overlap(:))
            ordering      = sort([p(1),q(1)]); % we merge to the first one
            merge_remove = ordering(2); 
            merge_keep   = ordering(1);

           
            [active, activity, phi, region,...
             timeSeries, cellIndex] =...
             merge(merge_keep, merge_remove, phi, c0, timeSeries,...
                   video_reshaped, t_len,...
                   mods, se_narrowband,...
                   active, activity,...
                   region, cellIndex); 

            p(1)               = [];
            q(1)               = [];
            p(p==merge_remove) = merge_keep;
            q(q==merge_remove) = merge_keep;
            p(p>merge_remove)  = p(p>merge_remove) - 1;
            q(q>merge_remove)  = q(q>merge_remove) - 1;
            unique_pq          = unique([p, q], 'rows');
            if ~isempty(unique_pq)
                p                  = unique_pq(:,1);
                q                  = unique_pq(:,2);
            end

        else
            p(1) = [];
            q(1) = [];
        end
    end    
end   


%%% Remove any that are too small or too big
nnzs                      = squeeze(sum(sum(phi<0,1),2));
outsized                  = or(nnzs<minimumSize, nnzs>maximumSize);
phi(:,:,outsized)         = [];
timeSeries(outsized,:)    = [];


disp(['Finished after ', num2str(it - 1), 'iterations']);
phi([1,end], :, :) = c0;
phi(:, [1,end], :) = c0;
cellMasks        = phi<0;    
 

function [] = updatePlot(ii, video, phi, radius, nrows, ncols, corrIm, meanIm,...
                         onlyThisMask, nhbd, region, diracPhi,...
                         similarityVelocity, regularisationVelocity,...
                         mu, mods, video_reshaped, t_len, timeSeries,...
                         it, activity, delta, noChangeNum)
     
    t = 1:size(video,3);
    [x,y] = find(phi(:,:,ii)<0);
    d = 2*radius;
    l_x = max(1, min(x) - d);
    u_x = min(size(video,1), max(x) + d);
    l_y = max(1, min(y) - d);
    u_y = min(size(video,2), max(y) + d);
     
    %%%% Mean image plot
    h(1) = subplot(nrows, ncols, 1);
    imagesc(corrIm(l_x:u_x,l_y:u_y));
    title('Corr image')
 
    %%%% Mean image and initial cell interior
    h(2) = subplot(nrows, ncols, 2);
    imagesc(meanIm(l_x:u_x,l_y:u_y));
    title('Mean image')
 
    %%%% Mean image and current cell interior
    h(3) = subplot(nrows, ncols, 3);   
    imshowpair(corrIm(l_x:u_x,l_y:u_y), phi(l_x:u_x,l_y:u_y,ii)<0)
    title(['Current interior ',num2str(ii)]);
 
    %%%% Mean image and current cell interior
    h(4) = subplot(nrows, ncols, 4);   
    imshowpair(corrIm(l_x:u_x,l_y:u_y), onlyThisMask(l_x:u_x,l_y:u_y))
    title(['Current solo interior ',num2str(ii)]);
 
     
    %%%% Mean image and narrowband around cell
    h(5) = subplot(nrows,ncols,5);   
    imshowpair(corrIm(l_x:u_x,l_y:u_y), nhbd(l_x:u_x,l_y:u_y))
    title('Current narrowband')
 
    %%%% Phi in 3D
    subplot(nrows,ncols,6);
    mesh(-phi(l_x:u_x,l_y:u_y,ii));   % for a better view, the LSF is displayed upside down
    hold on;  contour(phi(l_x:u_x,l_y:u_y,ii), [0,0], 'r','LineWidth',2);
    title('Current level set function');
 
    %%%% Other cell interiors
    if size(phi,3)>1
        h(6) = subplot(nrows, ncols, 7);   
        imshowpair(corrIm(l_x:u_x,l_y:u_y), any(phi(l_x:u_x,l_y:u_y,1:end~=ii)<0,3))
        title('Other cells (inactive and active)');
    end
     
    %%%% Cell similarity velocity    
    h(7) = subplot(nrows,ncols,8);
    V    = region(l_x:u_x,l_y:u_y,ii).*...
           diracPhi(l_x:u_x,l_y:u_y).*similarityVelocity(l_x:u_x,l_y:u_y);
    imagesc(V);
    title('-dirac * sim vel');
    colorbar;
    cmap           = cmocean('balance');
    u_lim = max(max(abs(V(:))),eps);
    set(gca,'CLim',[-u_lim u_lim])
    colormap(cmap);
    title('Cell based velocity');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
 
    %%%% Regularisation velocity
    h(8) = subplot(nrows,ncols,9);
    V = mu*region(l_x:u_x,l_y:u_y,ii).*...
           regularisationVelocity(l_x:u_x,l_y:u_y);
    imagesc(V);
    u_lim = max(max(abs(V(:))),eps);
    set(gca,'CLim',[-u_lim u_lim])
    colorbar;
    title('mu * Regularisation velocity');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
 
    %%%% Combined velocity
    h(9) = subplot(nrows, ncols, 10);
    V = region(l_x:u_x,l_y:u_y,ii).*...
                                    (mu * regularisationVelocity(l_x:u_x,l_y:u_y) - ...
                                     diracPhi(l_x:u_x,l_y:u_y).*...
                                     (similarityVelocity(l_x:u_x,l_y:u_y))); 
    imagesc(V);
    u_lim = max(max(abs(V(:))),eps);
    set(gca,'CLim',[-u_lim u_lim])
    colorbar;
    title('Combined velocity');
     set(gca,'xtick',[])
    set(gca,'ytick',[])    
 
    %%%% Time series plot: time series of cell interior and narrowband
    subplot(nrows,ncols, 2*ncols + (1:ncols));
    n_ts    = extractTimeSeries(region(:,:,ii), mods, video_reshaped, t_len);
    ts_plot = plotProperHeight( [timeSeries(ii,:)./norm(timeSeries(ii,:));...
                                n_ts./norm(n_ts)] );
    plot(t, ts_plot')
    hold on; 
    xlabel( 'Time(s)' );
    ylabel( 'DF/F' );
    legend('Interior','Narrowband');
    title(['Region time series, it = ',num2str(it), ' cell num = ',...
           num2str(ii)]);
    hold off;
 
    %%%% Activity plot: how many pixels joined and left contour
    subplot(nrows, ncols, ncols*3 + (1:ncols));
    plot(1:noChangeNum, activity(ii,:))
    hold on;
    line([0, noChangeNum],[delta, delta], 'LineStyle','- -','Color','r');
    xlabel( ['Previous ',num2str(noChangeNum),' iterations'] );
    ylabel( 'Number of pixels' );
    title(['No. of pixels on contour that changed, it = ',num2str(it)])
    hold off;
     
    linkaxes([h(1), h(2), h(3), h(4), h(5), h(6), h(7), h(8) ]);
    drawnow;
     
end
 
 
function [V] = cellSimilarityVelocity(video_dim, phi, n,...
                                    ii, narrowband,...
                                    mods,...
                                    video_reshaped,...
                                    t_len, timeSeries, overlapSaturation,...
                                    otherCells) 
    

    V                     = zeros(video_dim);
    a                     = squeeze(timeSeries(ii,:));
    narrowbandNoOverlap   = and(narrowband,~otherCells);
    narrowbandWithOverlap = and(narrowband, otherCells);
     
    % Calculate V(x) in pixels where there is no overlap
    if any(narrowbandNoOverlap(:))
        loc          = find(narrowbandNoOverlap);
        vid_loc      = bsxfun(@plus,loc, mods);
        raw_x        = single(video_reshaped(vid_loc(:)));
        raw_x_rs     = reshape(raw_x, length(loc), t_len);
        V_no_overlap = compare(raw_x_rs, n, a);
        V(narrowbandNoOverlap) = V_no_overlap;  
    end 
    
    % Calculate the velocity where there are overlapping cells
    if any(narrowbandWithOverlap(:))
         
        loc               = find(narrowbandWithOverlap);
        [loc_x,loc_y]     = ind2sub(video_dim,loc) ;
        vid_loc           = bsxfun(@plus,loc, mods);
        raw_x             = single(video_reshaped(vid_loc(:)));
        raw_x_rs          = reshape(raw_x, length(loc),t_len);
 
        for j = 1:length(loc) 
            J              = find(phi(loc_x(j),loc_y(j),:) < 0);
            J(J == ii)     = [];
            b              = sum(timeSeries(J,:),1);  % time series of cells that x is currently in
            x              = squeeze(raw_x_rs(j,:));     % time series of current pixel
 
            if overlapSaturation 
                a_plus_b = min(a+b,1);
            else
                a_plus_b = a + b;
            end
            
            %a_plus_b = a_plus_b/max(a_plus_b(:));
            
            V(loc(j))      = compare(x, b, a_plus_b);
            
        end
    end
end
 
 
function dist = compare(x, a_out, a_in)
  
    dist = zeros(size(x,1),1,'single');
   
   if strcmp(algorithmPars.metric, 'corr')
       dist = corr(x', a_in') - corr(x', a_out');
   elseif strcmp(algorithmPars.metric, 'euclid')
       dist = sum(a_out.^2) - sum(a_in.^2) + 2.* x * (a_in-a_out)';
       dist = dist./length(a_in);
   end
   
    
end
 
 
function V = distanceRegularisation(phi)
% This Matlab code implements an edge-based active contour model as an application of the Distance Regularized Level Set Evolution (DRLSE) formulation in Li et al's paper: C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation",  IEEE Trans. Image Processing, vol. 19 (12), pp.3243-3254, 2010.
% Author: Chunming Li, all rights reserved
% E-mail: lchunming@gmail.com   
%         li_chunming@hotmail.com 
% URL:  http://www.imagecomputing.org/~cmli/         
    [phi_x,phi_y] = gradient(phi);
 
    s   = sqrt(phi_x.^2 + phi_y.^2);
    a   = (s>=0) & (s<=1);
    b   = (s>1);
    ps  = a.*sin(2*pi*s)/(2*pi)+b.*(s-1);  % compute first order derivative of the double-well potential p2 in eqaution (16)
    dps = ((ps~=0).*ps+(ps==0))./((s~=0).*s+(s==0));  % compute d_p(s)=p'(s)/s in equation (10). As s-->0, we have d_p(s)-->1 according to equation (18)
    V   = div(dps.*phi_x - phi_x, dps.*phi_y - phi_y) + 4*del2(phi); 
    V([1,2,end-1,end], :)   = 0;
    V(:, [1,2,end-1,end])   = 0;
end
 
function f = div(nx,ny)
% This Matlab code implements an edge-based active contour model as an application of the Distance Regularized Level Set Evolution (DRLSE) formulation in Li et al's paper: C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation",  IEEE Trans. Image Processing, vol. 19 (12), pp.3243-3254, 2010.
% Author: Chunming Li, all rights reserved
% E-mail: lchunming@gmail.com   
%         li_chunming@hotmail.com 
% URL:  http://www.imagecomputing.org/~cmli/    
    [nxx,~]=gradient(nx);  
    [~,nyy]=gradient(ny);
    f=nxx+nyy;
end
 
function f = Dirac(x, sigma)
% This Matlab code implements an edge-based active contour model as an application of the Distance Regularized Level Set Evolution (DRLSE) formulation in Li et al's paper: C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation",  IEEE Trans. Image Processing, vol. 19 (12), pp.3243-3254, 2010.
% Author: Chunming Li, all rights reserved
% E-mail: lchunming@gmail.com   
%         li_chunming@hotmail.com 
% URL:  http://www.imagecomputing.org/~cmli/
f=(1/2/sigma)*(1+cos(pi*x/sigma));
b = (x<=sigma) & (x>=-sigma);
f = f.*b;
end
 
 
function g = NeumannBoundCond(f, video_dim)  
% This Matlab code implements an edge-based active contour model as an application of the Distance Regularized Level Set Evolution (DRLSE) formulation in Li et al's paper: C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation",  IEEE Trans. Image Processing, vol. 19 (12), pp.3243-3254, 2010.
% Author: Chunming Li, all rights reserved
% E-mail: lchunming@gmail.com   
%         li_chunming@hotmail.com 
% URL:  http://www.imagecomputing.org/~cmli/
    nrow = video_dim(1); ncol = video_dim(2);
    g = f;
    g([1 nrow],[1 ncol],:) = g([3 nrow-2],[3 ncol-2],:);  
    g([1 nrow],2:end-1,:) = g([3 nrow-2],2:end-1,:);          
    g(2:end-1,[1 ncol],:) = g(2:end-1,[3 ncol-2],:);  
end
 
function ts = extractTimeSeries(mask, mods, video_reshaped, t_len)
 
    mask([1,end], :) = 0;
    mask(:, [1,end]) = 0;
    loc              = find(mask);
    vid_loc          = bsxfun(@plus,loc, mods);
    raw_ts           = video_reshaped(vid_loc(:));
    raw_ts           = reshape(raw_ts, length(loc),t_len);
    ts               = mean(raw_ts,1, 'double'); %% much faster to calculate as double
    
end



function[active, activity, phi, region, timeSeries, cellIndex] =...
             merge(p, q, phi, c0, timeSeries, video_reshaped, t_len,...
                   mods, se_narrowband, active, activity,...
                   region,  cellIndex)
      
                
    % Update cell ii to be the union of cells ii and jj 
    phi(:,:,p)         = -c0 * any( phi(:,:,[p,q]) < 0, 3)+...
                          c0 *~any( phi(:,:,[p,q]) < 0, 3);
    timeSeries(p,:)    = extractTimeSeries( phi(:,:,p) < 0, mods,...
                                            video_reshaped, t_len);
    region(:,:,p)      = imdilate(phi(:,:,p)<0, se_narrowband);
     
    % Remove cell jj's information 
    phi(:,:,q)                  = [];
    timeSeries(q,:)             = [];
    active(q)                   = [];
    region(:,:,q)               = [];
    activity(cellIndex == q, :) = [];
    cellIndex                   = setdiff(cellIndex, q);
     
    if notifications_on
        disp(['ROIs ',num2str(p), ' and ', num2str(q),' merged.']);
    end
end
 
 
function[sd] = maskToSignedDistanceFunction(phi, c0)
 
sd        = double((phi > 0).*(bwdist(phi < 0)-0.5) - ...
                   (phi < 0).*(bwdist(phi > 0)-0.5));  
sd        = sd*c0;
sd(sd>c0) = c0;
sd(sd<-c0) = -c0;
 
end
 
end
