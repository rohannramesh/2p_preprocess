function [icaguidata] = rederive_icaguidata_with_selective_criterion_forSarah(MultiRunDirs,data,A,C,keep)
%%

% global icaguidata
sizY = data.sizY;

% to whittle down the ROIs
% this is the number of pixels in each poss ROI
% ind_all = throw_out_cells_too_small_too_large(keep,ROIvars);
% now lets throw out all things that don't meet that criterion
% [A,C,keep,ROIvars] = update_variables_for_ROIs_kept(ind_all,A,C,keep,ROIvars);
% lets throw out any ROIs in the ind_keep that are completely encircled in
% another ROI
% ind_all = eliminate_ROIs_with_overlap(A,keep,sizY);
% [A,C,keep,ROIvars] = update_variables_for_ROIs_kept(ind_all,A,C,keep,ROIvars);

% now check if there are masks from PCA/ICA and add them if they don't
% overlap
% first lets check for .signals file in each run
% [A,C,keep,ROIvars] = add_masks_from_PCA_ica(MultiRunDirs,sizY,A,C,keep,ROIvars);


% % do cross correlation of tc
% CC_tc = corr(ROIvars.tc');
% % this is the upper corner of the CC
% upperCC = triu(CC_tc+1,1);
% upperCC = upperCC(:);
% upperCC(upperCC == 0) = [];
% upperCC = upperCC - 1;


        




    

% this is what we are using to sort
% Sorting_metric = ROIvars.max_pr;
% SpatialW = ROIvars.rval_space;
% TemporalW = ROIvars.rval_time;
% now reset ind_keep
ind_keep = find(keep);
ind_nonkeep = find(keep == 0);
% these are the ROIs that the algorithm recommends keeping
icaguidata.suggested_ROI_keep = zeros(1,size(A,2));
icaguidata.ROIs_selected_from_orig = 1:length(keep);
% icaguidata.suggested_ROI_keep(ind_keep) = 1;
icaguidata.suggested_ROI_nonkeep = zeros(1,size(A,2));
% icaguidata.suggested_ROI_nonkeep(ind_nonkeep) = 1;


% simple fillers
icaguidata.countour_sd = 5;
icaguidata.fName = [MultiRunDirs.date_mouse '\temp.tif'];
icaguidata.masks = zeros(sizY(1),sizY(2));
icaguidata.masks_with_group = zeros(sizY(1),sizY(2));
icaguidata.masks_with_highcc = zeros(sizY(1),sizY(2));
icaguidata.ROI_list = zeros(1,size(A,2));
icaguidata.tc_cc = nan(size(A,2),size(A,2));
tc_tmp =  nan(size(C,2),size(A,2)); % for cc in the future
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
[MM,rrr] = sort(ind_keep,'ascend');
sorted_ind_keep = ind_keep(rrr);
for ii = 1:length(sorted_ind_keep)
    icaguidata.ica(indrr).filter = full(reshape(A(:,sorted_ind_keep(ii)),sizY(1),sizY(2)));
    icaguidata.ica(indrr).temporal_weights = C(sorted_ind_keep(ii),:);
    icaguidata.ica(indrr).trace = C(sorted_ind_keep(ii),:);
    % all masks
    % lets get centroid
    centroid_new = get_centroid(icaguidata.ica(indrr).filter);
    icaguidata.ica(indrr).centroid = centroid_new;
    icaguidata.ica(indrr).roiNum = [];
    [icaguidata.ica(indrr).idx(:,1) icaguidata.ica(indrr).idx(:,2)] = find(single(icaguidata.ica(indrr).filter>0));
    % lets also set this value to 1 in the ROI list variable
    icaguidata.ROI_list(indrr) = 1;
    icaguidata.suggested_ROI_keep(indrr) = 1;
    icaguidata.Sorting_metric(indrr) = 0;
    icaguidata.Weights.Spatial(indrr) = 0;
    icaguidata.Weights.Temporal(indrr) = 0;
%     icaguidata.Weights.Spatial(indrr) = SpatialW(sorted_ind_keep(ii));
%     icaguidata.Weights.Temporal(indrr) = TemporalW(sorted_ind_keep(ii));
    tc_tmp(:,indrr) = detrend(C(sorted_ind_keep(ii),:));
%     icaguidata.tc_cc(indrr,1:length(sorted_ind_keep)) = CC_tc(sorted_ind_keep(ii),sorted_ind_keep);
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

% now get the cc
CC_tc = corr(tc_tmp);
icaguidata.tc_cc = CC_tc;
CC_THRESHOLD = 0.65;
for curr_ROI_tmp = 1:indrr-1
    % make an image of correlated masks
    if sum(icaguidata.tc_cc(curr_ROI_tmp,curr_ROI_tmp+1:end) > CC_THRESHOLD) > 0
        icaguidata.masks_with_highcc = icaguidata.masks_with_highcc + ...
            single(icaguidata.ica(curr_ROI_tmp).filter>0);
    end
end
    

% now non classifed cells that should look at 
% if have a pr > 0.2 then highlight with a mask - for weeding out
thresh_inclusion_of_nonkeep = 1; % used to be 0.1 but not keeping any of these for now
[MM,rrr] = sort(ind_nonkeep,'ascend');
sorted_ind_nonkeep = ind_nonkeep(rrr);
for ii = 1:length(sorted_ind_nonkeep)
    icaguidata.ica(indrr).filter = full(reshape(A(:,sorted_ind_nonkeep(ii)),sizY(1),sizY(2)));
    icaguidata.ica(indrr).temporal_weights = C(sorted_ind_nonkeep(ii),:);
    icaguidata.ica(indrr).trace = C(sorted_ind_nonkeep(ii),:);
    % lets get centroid
    centroid_new = get_centroid(icaguidata.ica(indrr).filter);
    icaguidata.ica(indrr).centroid = centroid_new;
    icaguidata.ica(indrr).roiNum = [];
    [icaguidata.ica(indrr).idx(:,1) icaguidata.ica(indrr).idx(:,2)] = find(single(icaguidata.ica(indrr).filter>0));
    % for list purposes
    icaguidata.suggested_ROI_nonkeep(indrr) = 1;
    icaguidata.Sorting_metric(indrr) = 2;
    icaguidata.Weights.Spatial(indrr) = 2;
    icaguidata.Weights.Temporal(indrr) = 2;
%     icaguidata.Sorting_metric(indrr) = MM(ii);
%     icaguidata.Weights.Spatial(indrr) = SpatialW(sorted_ind_nonkeep(ii));
%     icaguidata.Weights.Temporal(indrr) = TemporalW(sorted_ind_nonkeep(ii));
%     icaguidata.tc_cc(indrr,1:length(sorted_ind_nonkeep)) = CC_tc(sorted_ind_nonkeep(ii),sorted_ind_nonkeep);
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
if ~isempty(ind_toss)
    ind_keep(ismember(ind_keep,ind_toss(:,1))) = [];
end
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

function [A,C,keep,ROIvars] = add_masks_from_PCA_ica(MultiRunDirs,sizY,A,C,keep,ROIvars);

% check for .signals file
have_previous_signals_file = 0;
for i = 1:length(MultiRunDirs.runs)
    if exist([MultiRunDirs.runs{i}.path '\' MultiRunDirs.runs{i}.sbx_name '.signals']) && ...
        have_previous_signals_file == 0
        have_previous_signals_file = 1;
        display('adding PCA/ICA masks')
        run_to_keep = i;
    end
end
if have_previous_signals_file == 1
    % load in signals file
    TMP = load([MultiRunDirs.runs{run_to_keep}.path '\' MultiRunDirs.runs{run_to_keep}.sbx_name '.signals'],'-mat');
    % make sparse matrix with all of masks
    PCA_ICAVar.Masks_PCA_ICA = sparse(prod(sizY(1:2)),length(TMP.cellsort));
    PCA_ICAVar.tc = nan(length(TMP.cellsort),length(TMP.cellsort(1).ica_trace));
    PCA_ICAVar.max_pr = nan(length(TMP.cellsort),1); 
    PCA_ICAVar.temporal_weights = nan(length(TMP.cellsort),1); 
    PCA_ICAVar.rval_space = nan(length(TMP.cellsort),1); 
    PCA_ICAVar.rval_time = nan(length(TMP.cellsort),1); 
    PCA_ICAVar.sizeA = nan(length(TMP.cellsort),1); 
    for curr_ROI = 1:length(TMP.cellsort)
        curr_mask = TMP.cellsort(curr_ROI).mask;
        % removed edges and ds
        curr_mask = sbxRemoveEdgesRR(curr_mask);
        curr_mask = bin_mov_xyt(curr_mask,2,1,0);
        % put in mask and tc
        PCA_ICAVar.Masks_PCA_ICA(:,curr_ROI) = sparse(reshape(curr_mask,prod(sizY(1:2)),[]));
        PCA_ICAVar.tc(curr_ROI,:) = TMP.cellsort(curr_ROI).ica_trace;
        PCA_ICAVar.sizeA(curr_ROI,1) = nnz(double(curr_mask>0));
    end
    % now lets test how many of the masks have large overlap
    options.dist_maxthr = 0.15;
    options.dist_exp = 1; 
    options.dist_thr = 0.5; 
    options.dist_overlap_thr = 0.8; 
    M1 = PCA_ICAVar.Masks_PCA_ICA > 0;
    M2 = A > 0;
    K1 = size(M1,2)-1; % for npil
    K2 = size(M2,2);
    D = zeros(K1,K2);
    for i = 1:K1
        for j = 1:K2

            overlap = nnz(M1(:,i) & M2(:,j));
            totalarea = nnz(M1(:,i)|M2(:,j));
            smallestROI = min(nnz(M1(:,i)),nnz(M2(:,j)));

            D(i,j) = 1 - (overlap/totalarea)^options.dist_exp;

            if overlap >= options.dist_overlap_thr*smallestROI
                D(i,j) = 0;
            end        

        end
    end
    D(D>options.dist_thr) = Inf;
    R = Hungarian(D);
    [match_1,match_2] = find(R);
%     match_1 = [];
    matched_ROIs = [match_1,match_2];
    % these are the ROI from PCA/ICA that aren't identified
    nonmatched_1 = setdiff(1:K1,match_1);
    nonmatched_2 = setdiff(1:K2,match_2);
    % set those that are matched via the algorithm to a 1 in keep
    keep(match_2) = 1;
    % now lets add them to the list - AT THE FRONT
    A = [PCA_ICAVar.Masks_PCA_ICA(:,nonmatched_1) A];
    C = [nan(length(nonmatched_1),size(C,2)); C];
    keep = [ones(length(nonmatched_1),1); keep];
    fields_ROIvars = fieldnames(ROIvars);
    for i = 1:length(fields_ROIvars)
        ROIvars.(fields_ROIvars{i}) = [PCA_ICAVar.(fields_ROIvars{i})(nonmatched_1,1:size(ROIvars.(fields_ROIvars{i}),2)); ...
            ROIvars.(fields_ROIvars{i})];
    end
else
    % just so don't get error
    RRR.A = A;
    RRR.C = C;
    RRR.keep = keep;
    RRR.ROIvars = ROIvars;
    PCA_ICAVar = [];
end
   
        
    

    
    