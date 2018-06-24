function [icaguidata] = convert_nmf_to_icaguidata(MultiRunDirs,data,A,C,keep,ROIvars)
%%

% global icaguidata
sizY = data.sizY;

ind_all = 1:length(keep);
% % to whittle down the ROIs
% % this is the number of pixels in each poss ROI
% ind_all = throw_out_cells_too_small_too_large(keep,ROIvars);
% % now lets throw out all things that don't meet that criterion
% [A,C,keep,ROIvars] = update_variables_for_ROIs_kept(ind_all,A,C,keep,ROIvars);
% % lets throw out any ROIs in the ind_keep that are completely encircled in
% % another ROI
% ind_all = eliminate_ROIs_with_overlap(A,keep,sizY);
% [A,C,keep,ROIvars] = update_variables_for_ROIs_kept(ind_all,A,C,keep,ROIvars);



    

% this is what we are using to sort
Sorting_metric = ROIvars.max_pr;
SpatialW = ROIvars.rval_space;
TemporalW = ROIvars.rval_time;
% now reset ind_keep
ind_keep = find(keep);
ind_nonkeep = find(keep == 0);
% these are the ROIs that the algorithm recommends keeping
icaguidata.suggested_ROI_keep = zeros(1,size(A,2));
icaguidata.ROIs_selected_from_orig = ind_all;
% icaguidata.suggested_ROI_keep(ind_keep) = 1;
icaguidata.suggested_ROI_nonkeep = zeros(1,size(A,2));
% icaguidata.suggested_ROI_nonkeep(ind_nonkeep) = 1;


% simple fillers
icaguidata.countour_sd = 5;
icaguidata.fName = [MultiRunDirs.date_mouse '\temp.tif'];
icaguidata.masks = zeros(sizY(1),sizY(2));
icaguidata.masks_with_group = zeros(sizY(1),sizY(2));
icaguidata.ROI_list = zeros(1,size(A,2));
% let make temporary mean mov 
fields_data = fieldnames(data);
if sum(ismember(fields_data,'movm'))
    icaguidata.movm = data.movm;
    icaguidata.movcorr = data.CrossCorrMeanImg;
else
    small_mov = data.Y(:,:,1:2000);
    icaguidata.movm = mean(small_mov,3);
    icaguidata.movcorr = mean(small_mov,3);
end

% preallocate SparseMasks
indrr = 1;
% first lets do sorted cells to keep
[MM,rrr] = sort(Sorting_metric(ind_keep),'descend');
sorted_ind_keep = ind_keep(rrr);
for ii = 1:length(sorted_ind_keep)
    icaguidata.ica(indrr).filter = full(reshape(A(:,sorted_ind_keep(ii)),sizY(1),sizY(2)));
    icaguidata.ica(indrr).temporal_weights = C(sorted_ind_keep(ii),:);
    icaguidata.ica(indrr).trace = ROIvars.tc(sorted_ind_keep(ii),:);
    % all masks
    % lets get centroid
    centroid_new = get_centroid(icaguidata.ica(indrr).filter);
    icaguidata.ica(indrr).centroid = centroid_new;
    icaguidata.ica(indrr).roiNum = [];
    [icaguidata.ica(indrr).idx(:,1) icaguidata.ica(indrr).idx(:,2)] = find(single(icaguidata.ica(indrr).filter>0));
    % lets also set this value to 1 in the ROI list variable
    icaguidata.ROI_list(indrr) = 1;
    icaguidata.suggested_ROI_keep(indrr) = 1;
    icaguidata.Sorting_metric(indrr) = MM(ii);
    icaguidata.Weights.Spatial(indrr) = SpatialW(sorted_ind_keep(ii));
    icaguidata.Weights.Temporal(indrr) = TemporalW(sorted_ind_keep(ii));
    % update the ROI masks image
        icaguidata.masks = icaguidata.masks + ...
            indrr*single(icaguidata.ica(indrr).filter>0);
        % masks with group number
        icaguidata.masks_with_group = icaguidata.masks_with_group + ...
            single(icaguidata.ica(indrr).filter>0);
        % can't have over 5 groups
        icaguidata.masks_with_group(icaguidata.masks_with_group > 6) = 6;
    indrr = indrr + 1;
end
% now non classifed cells that should look at 
% if have a pr > 0.2 then highlight with a mask - for weeding out
thresh_inclusion_of_nonkeep = 0.1;
[MM,rrr] = sort(Sorting_metric(ind_nonkeep),'descend');
sorted_ind_nonkeep = ind_nonkeep(rrr);
for ii = 1:length(sorted_ind_nonkeep)
    icaguidata.ica(indrr).filter = full(reshape(A(:,sorted_ind_nonkeep(ii)),sizY(1),sizY(2)));
    icaguidata.ica(indrr).temporal_weights = C(sorted_ind_nonkeep(ii),:);
    icaguidata.ica(indrr).trace = ROIvars.tc(sorted_ind_nonkeep(ii),:);
    % lets get centroid
    centroid_new = get_centroid(icaguidata.ica(indrr).filter);
    icaguidata.ica(indrr).centroid = centroid_new;
    icaguidata.ica(indrr).roiNum = [];
    [icaguidata.ica(indrr).idx(:,1) icaguidata.ica(indrr).idx(:,2)] = find(single(icaguidata.ica(indrr).filter>0));
    % for list purposes
    icaguidata.suggested_ROI_nonkeep(indrr) = 1;
    icaguidata.Sorting_metric(indrr) = MM(ii);
    icaguidata.Weights.Spatial(indrr) = SpatialW(sorted_ind_nonkeep(ii));
    icaguidata.Weights.Temporal(indrr) = TemporalW(sorted_ind_nonkeep(ii));
    if icaguidata.Sorting_metric(indrr) > thresh_inclusion_of_nonkeep
        group_showing = 3;
        icaguidata.ROI_list(indrr) = group_showing;
        % update the ROI masks image
        icaguidata.masks = icaguidata.masks + ...
            indrr*single(icaguidata.ica(indrr).filter>0);
        % masks with group number
        icaguidata.masks_with_group = icaguidata.masks_with_group + ...
            group_showing.*single(icaguidata.ica(indrr).filter>0);
        % can't have over 5 groups
        icaguidata.masks_with_group(icaguidata.masks_with_group > 6) = 6;
    end
    indrr = indrr + 1;
end   


function [A,C,keep,ROIvars] = update_variables_for_ROIs_kept(ind_keep,A,C,keep,ROIvars);

A = A(:,ind_keep);
C = C(ind_keep,:);
keep = keep(ind_keep);
fields_ROIvars = fieldnames(ROIvars);
for i = 1:length(fields_ROIvars)
    ROIvars.(fields_ROIvars{i}) = ROIvars.(fields_ROIvars{i})(ind_keep,:);
end

function ind_all = eliminate_ROIs_with_overlap(A,keep,sizY)

ind_keep = find(keep);
tmp_keep_masks = zeros(sizY(1),sizY(2),length(ind_keep));
for i = 1:length(ind_keep)
    tmp_keep_masks(:,:,i) = full(reshape(A(:,ind_keep(i)),sizY(1),sizY(2))>0);
end
% now ask for each ROI is the entirety of the ROI in another ROI

ind_toss = [];
for i = 1:length(ind_keep)
    curr_mask1 = tmp_keep_masks(:,:,i);
    pix_mask1 = find(curr_mask1);
    other_ROIs = setdiff(1:length(ind_keep),i);
    other_ROIs(other_ROIs < i) = [];
    for j = other_ROIs;
        curr_mask2 = tmp_keep_masks(:,:,j);
        pix_mask2 = find(curr_mask2);
        % this is the % overlap
        if ~isempty(intersect(pix_mask1,pix_mask2))
            Percover(1) = length(intersect(pix_mask1,pix_mask2))./length(pix_mask1);
            Percover(2) = length(intersect(pix_mask2,pix_mask1))./length(pix_mask2);
            % if overlap greater than 90% then toss that cell
            thresh = 0.9;
            if sum(Percover > thresh) == 1
                % which ROI has the highest overlap
                tmpR = find(Percover > thresh);
                if tmpR == 1
                    ind_toss = [ind_toss; ind_keep(i) ind_keep(j)];
                elseif tmpR == 2
                    ind_toss = [ind_toss; ind_keep(j) ind_keep(i)];
                end  
             % if both cells have overlap > thresh keep bigger ROI
            elseif sum(Percover > thresh) == 2                
                nnz_pix = [length(pix_mask1) length(pix_mask2)];
                tmpR = find(nnz_pix == max(nnz_pix))
                if tmpR == 1
                    ind_toss = [ind_toss; ind_keep(i)  ind_keep(j)];
                elseif tmpR == 2
                    ind_toss = [ind_toss; ind_keep(j)  ind_keep(i)];
                end
            end
        else
            Percover = [0 0];
        end
    end
end
% now remove the ind_toss first column from keep
ind_keep(ismember(ind_keep,ind_toss(:,1))) = [];
ind_nonkeep = find(keep == 0);
ind_all = sort([ind_keep; ind_nonkeep]);

function ind_all = throw_out_cells_too_small_too_large(keep,ROIvars);

% now put things in proper order for the icaguidata
% first number is the ones to keep
ind_keep = find(keep);
ind_nonkeep = find(keep == 0);

min_ROI_size = 7; % in pixels
max_ROI_size = 500;
ind_too_small = find(ROIvars.sizeA <= min_ROI_size | ROIvars.sizeA >= max_ROI_size);
% remove any ROI that are less than 5 pixels
ind_keep(ismember(ind_keep,ind_too_small)) = [];
ind_nonkeep(ismember(ind_nonkeep,ind_too_small)) = [];
ind_all = sort([ind_keep; ind_nonkeep]);
