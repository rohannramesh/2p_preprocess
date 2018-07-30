function all_task_variables = get_task_variables_for_GLM(Dirs, variables_to_consider, ds_val)
% for all variables in variables to consider create a vector that
% represents that task variable for the future model and build a matrix of
% these variables

if nargin < 3
    ds_val = 1;
end


% iterate through each run and get all of the variables we care about
nRuns = length(Dirs.runs);
% to concatenate across runs
all_task_variables = [];
for curr_run = 1:nRuns
    % this is the directory information for a specific run
    curr_dir = Dirs.runs{curr_run};
    % load in the frame_2p_metadata and the TrialVariable
    load(strrep(curr_dir.sbx,'.sbx','.f2p'),'-mat');
    load(strrep(curr_dir.sbx,'.sbx','.TrialVar'),'-mat');
    

    % iterate through each variable and load that vector
    all_task_variables_curr_run = [];
    for i = 1:length(variables_to_consider)
        curr_task_variable = variables_to_consider{i};
        
        % get the appropriate task variable
        output = process_task_variables(frame_2p_metadata,Stim_Master,Column_Order,curr_task_variable,curr_dir,ds_val);
        
        % make sure it is the proper length
        output = output(1:floor(length(frame_2p_metadata.ensuretrigsperframe)/ds_val));
        
        % make sure it is proper orientation for concatenations
        if size(output,1) > 1
            output = output';
        end
        
        % concatenate
        all_task_variables_curr_run = [all_task_variables_curr_run; output];  
    end
    % to concatenate across runs
    all_task_variables = [all_task_variables all_task_variables_curr_run];
end

end


function output = process_task_variables(frame_2p_metadata,...
    Stim_Master,Column_Order,curr_task_variable,curr_dir,ds_val)

% if lick bout onset then only take the first lick in a lickbout
if strcmp('lick bout onset',curr_task_variable) | strcmp('all other licking',curr_task_variable)
    licking = frame_2p_metadata.lickingtrigsperframe;
    output = process_licking(licking,curr_task_variable,frame_2p_metadata);
    if ds_val ~= 1
        output = PPPack.hf.downsample_vector(output,ds_val,'max');
    end
elseif strcmp('quinine',curr_task_variable)
    output = frame_2p_metadata.quininetrigsperframe;
    if ds_val ~= 1
        output = PPPack.hf.downsample_vector(output,ds_val,'sum');
    end
elseif strcmp('ensure',curr_task_variable)
    output = frame_2p_metadata.ensuretrigsperframe;
    if ds_val ~= 1
        output = PPPack.hf.downsample_vector(output,ds_val,'sum');
    end
elseif strcmp('xyshift',curr_task_variable)
    % load in the .align file
    TMP = load(strrep(curr_dir.sbx,'.sbx','.align'),'-mat');
    xshift = TMP.StackReg(:,4);
    yshift = TMP.StackReg(:,3);
    output = sqrt(xshift.^2 + yshift.^2); % this is the brain movement variable
    output = diff([output(1); output]); % want to look at changes so do diff
    if ds_val ~= 1
        output = PPPack.hf.downsample_vector(output,ds_val,'mean');
    end
elseif ~isempty(strfind(curr_task_variable,'onsets')) || ...
        ~isempty(strfind(curr_task_variable,'entire')) || ...
        ~isempty(strfind(curr_task_variable,'offsets')) % for all visual stimulus stuff - giant function below
    output = build_visual_task_variable(frame_2p_metadata,Stim_Master,Column_Order,curr_task_variable);
    if ds_val ~= 1
        output = PPPack.hf.downsample_vector(output,ds_val,'max');
    end
elseif strcmp('running',curr_task_variable)   
    [mouse,curr_date,run] = PPPack.hf.get_mouse_day_run_info_from_dirs(curr_dir);
    output = PPPack.hf.sbxSpeed(mouse,curr_date,run);
    if ds_val ~= 1
        output = PPPack.hf.downsample_vector(output,ds_val,'mean');
    end
end

end


function new_lick_vec = process_licking(licking,curr_task_variable,frame_2p_metadata)
    % define time without a lick to define between bout onsets
    min_time.sec = 1;
    min_time.fr = round(min_time.sec*frame_2p_metadata.meta.framerate_2p);
    ts_each_lick = find(licking);
    % find those timestamps that are more than the min apart
    diff_ts = diff([-50 ts_each_lick]);
    % big gap
    ts_idx_w_min_time = find(diff_ts > min_time.fr);
    % now to make the vector
    lick_onsets_vec = zeros(size(licking));
    lick_onsets_vec(ts_each_lick(ts_idx_w_min_time)) = 1;
    % for all other licks
    all_other_licks_vec = zeros(size(licking));
    all_other_licks_vec(setdiff(ts_each_lick,ts_each_lick(ts_idx_w_min_time))) = 1;
    if strcmp('lick bout onset',curr_task_variable)
        new_lick_vec = lick_onsets_vec;
    elseif strcmp('all other licking',curr_task_variable)
        new_lick_vec = all_other_licks_vec;
    end
end

function iOnsets_inds = get_iOnsets_stimtype(conditions_name,dataType,conditions_to_use)
% get the indices of onsets of a particular stimulus type
if nargin < 3
    conditions_to_use = {'CSp','CSm_cond','CSn_cond'}; 
end

% these are the condition numbers
unique_conditions = unique(conditions_name,'stable');
% lets identify which conditions want to keep
iOnsets_inds = zeros(size(conditions_name));
if strfind(dataType,'FC')
    iOnsets_inds = [iOnsets_inds+ ...
    (cellfun(@isempty,strfind(conditions_name,conditions_to_use{1}))==0)];
elseif strfind(dataType,'QC')
    iOnsets_inds = [iOnsets_inds+...
    (cellfun(@isempty,strfind(conditions_name,conditions_to_use{2}))==0)];
elseif strfind(dataType,'NC')
    iOnsets_inds = [iOnsets_inds+ ...
    (cellfun(@isempty,strfind(conditions_name,conditions_to_use{3}))==0)];
end
iOnsets_inds = iOnsets_inds>0;
end

function iOnsets_inds = get_iOnsets_behaviortype(conditions_name,trialerror,dataType)
% to get indices of onsets of a particular behavior type, i.e. did they
% perform the trial correctly or incorrectly
% these are the condition numbers
unique_conditions = unique(conditions_name,'stable');
% lets identify which conditions want to keep
curr_te = trialerror;
iOnsets_inds = zeros(size(curr_te));
if strfind(dataType,'incorrect')
    iOnsets_inds = [iOnsets_inds+ mod(curr_te,2)==1];
elseif strfind(dataType,'correct')
    iOnsets_inds = [iOnsets_inds+ mod(curr_te,2)==0];
end
iOnsets_inds = iOnsets_inds>0;
end

function curr_conditions_name = build_conditions_name_variable(frame_2p_metadata,Stim_Master,Column_Order)
% build a cell array so that you can identify exactly
% lets make conditions name variable
curr_conditions = Stim_Master(:,ismember(Column_Order,'Condition')); % seventh column corresponds to stimulus condition number
curr_conditions_name = cell(size(curr_conditions));
if isfield(frame_2p_metadata.meta,'condition_number')
    for curr_cond_ind = 1:length(frame_2p_metadata.meta.condition_number)
        indtmp = find(curr_conditions == frame_2p_metadata.meta.condition_number(curr_cond_ind));
        for i = indtmp'
            curr_conditions_name{i} = frame_2p_metadata.meta.condition_name{curr_cond_ind};
        end
    end
end
end

function output = build_visual_task_variable(frame_2p_metadata,Stim_Master,Column_Order,curr_task_variable)
    % build the conditions name variable for indexing of stimulus type
    conditions_name = build_conditions_name_variable(frame_2p_metadata, Stim_Master, Column_Order);
    % get the stimulus offsets
    logic    = (frame_2p_metadata.contrast ~= 0);
    difflogic = [0 diff(logic)];
    iOnsets  = find(difflogic == 1);
    iOffsets = find(difflogic == -1);
    iOffsets(iOffsets<iOnsets(1)) = [];
    iOffsets = iOffsets - 1;
    if frame_2p_metadata.contrast(end) ~= 0
        iOffsets(end+1) = ...
            find(isnan(frame_2p_metadata.contrast)==0, 1, 'last');
    end
    % lets do for stimulus onsets
    if ~isempty(strfind(curr_task_variable,'onsets'))
        % if FC QC or NC subdivide by the cue type
        if ~isempty(strfind(curr_task_variable,'FC')) || ...
            ~isempty(strfind(curr_task_variable,'QC')) || ...
            ~isempty(strfind(curr_task_variable,'NC')) 
                iOnsets_inds = get_iOnsets_stimtype(conditions_name,curr_task_variable);
                output = zeros(size(frame_2p_metadata.quininetrigsperframe));
                tmpR = Stim_Master(:,ismember(Column_Order,'iOnsets'));
                tmpR = tmpR(find(iOnsets_inds));
                output(tmpR) = 1;
        else % if keeping all stimulus onsets
                output = zeros(size(frame_2p_metadata.quininetrigsperframe));
                tmpR = Stim_Master(:,ismember(Column_Order,'iOnsets'));
                output(tmpR) = 1;
        end
    elseif ~isempty(strfind(curr_task_variable,'offsets')) % for stimulus offsets
        % if FC QC or NC subdivide by the cue type
        if ~isempty(strfind(curr_task_variable,'FC')) || ...
            ~isempty(strfind(curr_task_variable,'QC')) || ...
            ~isempty(strfind(curr_task_variable,'NC')) 
                iOnsets_inds = get_iOnsets_stimtype(conditions_name,curr_task_variable);
                output = zeros(size(frame_2p_metadata.quininetrigsperframe));
                tmpR = iOffsets;
                tmpR = tmpR(find(iOnsets_inds));
                output(tmpR) = 1;
        else % if keeping all stimulus onsets
                output = zeros(size(frame_2p_metadata.quininetrigsperframe));
                tmpR = iOffsets;
                output(tmpR) = 1;
        end
    elseif ~isempty(strfind(curr_task_variable,'entire')) % for entire stimulus
        % if only want subportion of either FC QC or NC
        if ~isempty(strfind(curr_task_variable,'FC')) || ...
                ~isempty(strfind(curr_task_variable,'QC')) || ...
                ~isempty(strfind(curr_task_variable,'NC')) 
            % if sorting by behavioral performance
            if ~isempty(strfind(curr_task_variable,'incorrect')) || ~isempty(strfind(curr_task_variable,' correct')) 
                    iOnsets_inds = get_iOnsets_stimtype(conditions_name,curr_task_variable);
                    % correct onsets
                    if ~isempty(strfind(curr_task_variable,'incorrect')) 
                        iOnsets_inds_behav = get_iOnsets_behaviortype(conditions_name,...
                            Stim_Master(:,ismember(Column_Order,'TrialError')),'incorrect');
                    elseif ~isempty(strfind(curr_task_variable,' correct')) 
                        iOnsets_inds_behav = get_iOnsets_behaviortype(conditions_name,...
                            Stim_Master(:,ismember(Column_Order,'TrialError')),' correct');
                    end
                    output = zeros(size(frame_2p_metadata.quininetrigsperframe));
                    tmpR = Stim_Master(:,ismember(Column_Order,'iOnsets'));
                    tmpR2 = iOffsets;
                    % select only of stim type
                    tmpR = tmpR(intersect(find(iOnsets_inds),find(iOnsets_inds_behav)));
                    tmpR2 = tmpR2(intersect(find(iOnsets_inds),find(iOnsets_inds_behav)));
                    stim_entire = [];
                    for iii = 1:length(tmpR)
                        stim_entire = [stim_entire tmpR(iii):tmpR2(iii)];
                    end
                    output(stim_entire) = 1;           
            else % if not take all trials of a particular cue type
                    iOnsets_inds = get_iOnsets_stimtype(conditions_name,curr_task_variable);
                    output = zeros(size(frame_2p_metadata.quininetrigsperframe));
                    tmpR = Stim_Master(:,ismember(Column_Order,'iOnsets'));
                    tmpR2 = iOffsets;
                    % select only of stim type
                    tmpR = tmpR(find(iOnsets_inds));
                    tmpR2 = tmpR2(find(iOnsets_inds));
                    stim_entire = [];
                    for iii = 1:length(tmpR)
                        stim_entire = [stim_entire tmpR(iii):tmpR2(iii)];
                    end
                    output(stim_entire) = 1;
            end
        else % if taking all stimulus durations
            output = zeros(size(frame_2p_metadata.quininetrigsperframe));
            tmpR = Stim_Master(:,ismember(Column_Order,'iOnsets'));
            tmpR2 = iOffsets;
            stim_entire = [];
            for iii = 1:length(tmpR)
                stim_entire = [stim_entire tmpR(iii):tmpR2(iii)];
            end
            output(stim_entire) = 1;
        end
    end
end