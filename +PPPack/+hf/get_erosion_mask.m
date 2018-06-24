function curr_mask = get_erosion_mask(curr_mask, erosion_value)

curr_mask_norm = curr_mask./max(curr_mask(:));
baseline_threshold = 0;
nonzero_values = (curr_mask_norm(curr_mask_norm > 0));
[nonzero_values_sort sort_ind] = sort(nonzero_values,'ascend');


% turn erosion level into number of pixels
curr_threshold = floor((1-erosion_value)*length(nonzero_values_sort));
value_of_nonzero_values_to_toss = nonzero_values_sort(1:min([curr_threshold length(nonzero_values_sort)]));
curr_mask(ismember(curr_mask_norm,value_of_nonzero_values_to_toss)) = 0;


% additional step where ask if not connected to highest point then chuck
% image
% this is the largest point
try [idx,idy] = find(curr_mask == max(curr_mask(:)));
    idx = idx(1);idy = idy(1);
    % now see what are connected things
    bwc = bwconncomp(curr_mask>0,8);
    % make image of non connected things
    labmat = labelmatrix(bwc);
    % this is the value in the labelmatrix of the highest weight
    curr_labmat_highestpx = double(labmat(idx,idy));
    new_idx_to_toss = find(~ismember(labmat,[0 curr_labmat_highestpx]));
    curr_mask(new_idx_to_toss) = 0;
catch err
    warning('Struggle with highest point erosion')
end
end