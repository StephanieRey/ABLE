function [image] = crossCorr(video)

% Written by Stephanie Reynolds, 25/09/2017
% At each pixel this function outputs the average correlation of that
% pixel's time series and with the time series' of all the pixels in their 
% 8-connected neighbourhood.
%
%%%%%%%%%%%%%%%   INPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% video        MxNxT array
%
%%%%%%%%%%%%%%%   OUTPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%image         MxN array, the value of each pixel represents the average 
%              correlation for that pixel with its neighbours

dim             = [size(video,1), size(video,2)];
t_len           = size(video,3);
mods            = (0:1:(t_len - 1))*dim(1)*dim(2);
image           = zeros(dim);
video_reshaped  = reshape(video,[dim(1)*dim(2)*t_len,1,1]);

for ii = 1:dim(1)
    for jj = 1:dim(2)
            
            sub_loc      = [ii-1, jj-1;...
                            ii-1, jj;...
                            ii-1, jj+1;...
                            ii,   jj-1;...
                            ii,   jj+1;...
                            ii+1, jj-1;...
                            ii+1, jj;...
                            ii+1, jj+1];
            sub_loc      = get_locs_in_range(sub_loc, dim);
            sub_loc      = [ii, jj; sub_loc];
            ind_loc      = sub2ind(dim, sub_loc(:,1), sub_loc(:,2));
            vid_loc      = bsxfun(@plus, ind_loc, mods);
            raw_ts       = single(video_reshaped(vid_loc(:)));
            raw_ts_rs    = reshape(raw_ts, length(ind_loc), t_len);
            Y            = nanmean(raw_ts_rs([2:end],:),1);
            X            = raw_ts_rs(1,:);
            mu_Y         = mean(Y,2);
            sig_Y        = std(Y, [], 2);
            Y            = bsxfun(@minus, Y, mu_Y);
            Y            = bsxfun(@rdivide, Y, sig_Y);
            X            = (X - mean(X))/std(X);
            image(ii,jj) = Y*X';
            
    end
end

image = image./t_len;


function[sub_loc] = get_locs_in_range(sub_loc, dim)
    
    valid       = zeros(size(sub_loc));
    valid(:, 1) = and(sub_loc(:,1)>0, sub_loc(:,1) <= dim(1));
    valid(:, 2) = and(sub_loc(:,2)>0, sub_loc(:,2) <= dim(2));
    valid       = and(valid(:,1), valid(:,2));
    sub_loc     = sub_loc(valid, :);
    
end


end