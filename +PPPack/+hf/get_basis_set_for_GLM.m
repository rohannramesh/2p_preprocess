function [All_Basis_Sets,convolution_gauss,Beta_coeff_idx] = get_basis_set_for_GLM(Dirs, variables_to_consider, task_variables,...
    ds_val, prepost_time_windows_taskvar)

% gaussian to convolved for each basis set
convolution_gauss = gausswin(15); % std of about 1 sec if ds_val = 3 % used to be 15
convolution_gauss = convolution_gauss./sum(convolution_gauss);

% load in the frame_2p_metadata and the TrialVariable
load(strrep(Dirs.runs{1}.sbx,'.sbx','.f2p'),'-mat');
framerate = frame_2p_metadata.meta.framerate_2p;


% downsampling = 1;
nFrames_per_gauss = floor(ds_val*1.5);
if ds_val ~= 3
    warning('Make sure to check the gaussians in get_basis_set_for_GLM to see if overlap properly')
end

% iterate through each run and get all of the variables we care about
All_Basis_Sets = [];
Beta_coeff_idx.values = [];
for i = 1:length(variables_to_consider)
    % this is the current task variable and the vector associated with it
    curr_task_variable = variables_to_consider{i};
    vec = task_variables(i,:);
    
    % get the basis set for that variable
    Curr_Basis_Set = build_basis_set(vec,prepost_time_windows_taskvar(i,:),framerate,...
        ds_val,nFrames_per_gauss,curr_task_variable,convolution_gauss);
    
    % concatenate
    All_Basis_Sets = [All_Basis_Sets; Curr_Basis_Set'];
    Beta_coeff_idx.values = [Beta_coeff_idx.values; i*ones(size(Curr_Basis_Set',1),1)];
    Beta_coeff_idx.labels{i} = curr_task_variable;
    
end


end



function Curr_Basis_Set = build_basis_set(vec,prepost_time_extend,framerate,...
    downsampling,nFrames_per_gauss,curr_task_variable,w2);

% this is to deal with square wave pulses that last the duration of the stimulus
if ~isempty(strfind(curr_task_variable,' entire'))
    % reget onsets and offsets
    stim_onset = find(diff([0 vec]) == 1);
    stim_offset = find(diff([0 vec]) == -1);
    if length(stim_onset) > length(stim_offset)
        stim_onset(end) = [];
    end
    entire_stim = vec;
    % lets also extend stimulus by certain number of frames to
    % get decay - IGNORING INPUTS TO EXTEND BEFOREHAND BECAUSE DOESN"T MAKE
    % SENSE
    nExtend = round(prepost_time_extend(2)*framerate/downsampling); % in sec
    stim_offset = stim_offset + nExtend -1;
    for i = 1:length(stim_onset)
        entire_stim(stim_onset(i):stim_offset(i)) = 1;
    end  
    Curr_Basis_Set = [];
    circshift_vec = 0;
    while sum(entire_stim) > 0
        curr_behav_mat = zeros(size(vec));
        idx_with_val_of_one = find(diff([0 entire_stim]) == 1);
        curr_behav_mat(idx_with_val_of_one) = 1;
        curr_behav_mat = conv(curr_behav_mat,w2,'same');
        % have to throw out the onset now so the diff gets the
        % next frame
        for rr = 0:nFrames_per_gauss
            entire_stim(idx_with_val_of_one+rr) = 0;
        end
        Curr_Basis_Set = cat(2,Curr_Basis_Set,curr_behav_mat');
    end
elseif ~isempty(strfind(curr_task_variable,'xyshift')) || ~isempty(strfind(curr_task_variable,'running')) || ...
       ~isempty(strfind(curr_task_variable,'xshift')) || ~isempty(strfind(curr_task_variable,'yshift')) || ...
       ~isempty(strfind(curr_task_variable,'all other licking'))
    output = vec;
    Curr_Basis_Set = conv(output(:),w2,'same');
else % this is for any pulse extend some amount forward or backwards
    nPre = round(prepost_time_extend(1)*framerate/downsampling); % in sec  used to be 1 
    nPost = round(prepost_time_extend(2)*framerate/downsampling); % in sec used to be 2
    % put gaussian every this number of frames 
    % nFrames_per_gauss = 3;
    % this is vector of each fr pre and post target frame to use
    circshift_vec = [fliplr([0:-1*nFrames_per_gauss:-1*nPre]) ...
        ([nFrames_per_gauss:nFrames_per_gauss:nPost])];
    nTp_for_gaussian = length(circshift_vec);
    % Pre-allocate new basis set
    Curr_Basis_Set = zeros(length(vec),...
            nTp_for_gaussian);
    for i = 1:nTp_for_gaussian
        curr_behav_mat = zeros(size(vec));
        idx_with_val_of_one = find(vec);
        idx_with_val_of_one = idx_with_val_of_one+circshift_vec(i);
        idx_with_val_of_one(idx_with_val_of_one <= 0) = [];
        idx_with_val_of_one(idx_with_val_of_one > length(curr_behav_mat)) = [];
        curr_behav_mat(idx_with_val_of_one) = 1;
        curr_behav_mat = conv(curr_behav_mat,w2,'same');
        Curr_Basis_Set(:,i) = curr_behav_mat;
    end
end
end

