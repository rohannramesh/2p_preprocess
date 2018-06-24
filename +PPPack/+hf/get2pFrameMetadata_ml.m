function frame_2p_metadata...
	= get2pFrameMetadata_ml(stim, nidaq_data, n2pFrames,ch)

% get2pFrameMetadata takes the On the ball BHV script and NiDAQ-recorded timing
% data for a 2p movie and returns a struct with metadata for each frame in
% the movie.
%
% INPUTS:
%	stim		bhv file from MonkeyLogic (Asaad and Eskander, 2011)
%	nidaq_data	Timing data recorded via NiDAQ during 2p acquisition.
%               See below for channel numbers
%
%	n2pFrames	The recorded number of 2p frames. The number of frame
%				pulses that appear in the NiDAQ data are a few more than 
%				the recorded number of frames. Testing showed that these
%				extra pulses occur at the end of 2p acquisition such that
%				the first recorded pulse corresponds to the first recorded
%				TIFF frame.

%	
% OUTPUTS:
%	frame_2p_metadata	The structure that contains all the metadata for
%						the 2p frames.
%%




%% Define Conditions
if ~isempty(stim)
    Condition = stim.ConditionNumber; % Vector that lists which condition occurred each trial
    stim.frame.orientation = nan(size(Condition));
    stim.frame.sf          = nan(size(Condition));
    stim.frame.tf          = nan(size(Condition));
    stim.frame.location    = nan(size(Condition));
    stim.frame.condition   = Condition;
    stim.frame.trialerror  = stim.TrialError;
    stim.frame.block       = stim.BlockNumber;
    stim.frame.contrast    = nan(size(Condition));

    unique_Condition = unique(Condition);


    % Deal with an old formatting issue from older .bhv file types
    if ~isfield(stim,'RunType')
        RunType = 'Condition'; % options = Condition, Ori, Retinotopy
    elseif strcmp(stim.RunType,'Ori')
        RunType = 'Ori'; 
    elseif strcmp(stim.RunType,'Retinotopy')
        RunType = 'Retinotopy'; 
    elseif strcmp(stim.RunType,'Training')
        RunType = 'Condition'; 
        warning('Assuming Condition run');
    end



    % to parse information about the stimulus identity need a pointer to either
    % the Timing Scripts used by MonkeyLogic or to the TaskObjects file which
    % includes information about the movies presented to the mouse
    switch RunType
        case 'Condition'
            Timing_Scripts = stim.TimingFileByCond(unique_Condition);
        case 'Ori'
            Timing_Scripts = stim.TaskObject(unique_Condition,2); %for ori runs
        case 'Retinotopy'
            Timing_Scripts = stim.TaskObject(unique_Condition,2); %for retinotopy runs
        otherwise
            error('Not valid choice')
    end


    for count = 1:length(unique_Condition);
        curr_Timing_Script = Timing_Scripts(count);
        % Assumption: CSp = 0, CSn = 135; CSm = 270; Blank = 360 sf = 0.04, tf = 2, 
        % location = -1
        switch RunType
        case 'Condition'
            stim = PPPack.hf.parse_bhv_mov_info(stim,curr_Timing_Script,count);
        case 'Ori'
            stim = parse_orientation_tuning_movie(stim,Condition,unique_Condition,curr_Timing_Script,count)
        case 'Retinotopy'
            stim = parse_retinotopy_movies(stim,Condition,unique_Condition,curr_Timing_Script,count)
        end
    end
end

%%


multiWaitbar('Processing frame metadata', 0);

%% Sort Nidaq data
% The Nidaq data -- unbelievably -- is not saved in chronological order by
% the recording device. This is not an issue with regards to accuracy
% because timestamps are saved (correctly) for each data point. However, 
% We rely on chronological order below, so we have to sort the data here.
% This does not change the timing of each data point, it just changes the 
% order in which the data are saved in the arrays.
[nidaq_data.timestamps, time_idx] = sort(nidaq_data.timestamps);
for i = 1:size(nidaq_data.data, 2)
	nidaq_data.data(:, i) = nidaq_data.data(time_idx, i);
end

%% Extract running eyetracking and bodytracking events
% now for running
Nidaq_rate = 1./diff(nidaq_data.timestamps(1:2));
frame_running_onsets = get_onsets(nidaq_data,ch.running);


% for eye camera pulses
[frame_eyetrack_onsets,ind] = get_onsets(nidaq_data,ch.eyetrack);
frame_eyetrack_pulses = zeros(length(nidaq_data.data(:,ch.eyetrack)),1);
frame_eyetrack_pulses(ind) = 1;

% for licking
clear ind
frame_licking_onsets = get_onsets(nidaq_data,ch.licking);

% for ensure
clear ind
frame_ensure_onsets = get_onsets(nidaq_data,ch.ensure);

% for quinine
clear ind
frame_quinine_onsets = get_onsets(nidaq_data,ch.quinine);

% for body track
frame_bodytrack_onsets = get_onsets(nidaq_data,ch.bodytrack);


%% Extract times of stimulation frames (monitor framerate) from rec_data:

% Use simple thresholding to get frame onset times from the analog frame
% signal recorded by the NiDAQ system:
clear ind
frame_2p_onsets = get_onsets(nidaq_data,ch.twoP);


% if the number of f2p onsets is half the number of total movie frames
% this is bc the new 2p rig sends out only one pulse for bidirectional
% scanning
if abs(length(frame_2p_onsets)*2 - n2pFrames) <= 100 % need or bc sometimes ends in up state
    frame_2p_onsets = interp1(1:2:length(frame_2p_onsets)*2,frame_2p_onsets,0:length(frame_2p_onsets)*2-1,'linear','extrap');
    warning('Assuming on new 2p rig. Bc only send out one pulse for forward and back scan, artificially adding pulse for flyback.')
end
    
% if the number of 2p onsets doesn't match the number of frames give a
% warning
if numel(frame_2p_onsets) < n2pFrames
	w = warndlg(['Number of 2p onsets ' num2str(length(frame_2p_onsets)) ...
        ' number of frames in movie ' num2str(n2pFrames)]);
elseif numel(frame_2p_onsets) > n2pFrames
        pulse_extra=numel(frame_2p_onsets) - n2pFrames;
    w = warndlg(cat(2,'It looks like the Nidaq timing have ',num2str(pulse_extra), ' more pulses than movie frames.'));
end


if ~isempty(stim)
% Visstim pulses are +1V for all frames and +5V for non-blank frames
clear ind
thresh_stim = 0.5; % Get all frames by thresholding at 0.5 V
[frame_stim_onsets,ind] = get_onsets(nidaq_data,ch.OTB_visstim,thresh_stim);


% never allow a frame_stim_onset to occur within given time or frame window
min_time_between = 5;
diff_ind = diff([-1*min_time_between ;frame_stim_onsets]); % used to be median(diff(frame_stim_onsets))
ind_error = find(diff_ind < min_time_between);
frame_stim_onsets(ind_error) = [];

% get the frame stim offsets
ind2 = find(diff(nidaq_data.data(:,ch.OTB_visstim)>thresh_stim) == -1);
Nidaq_rate = 1./diff(nidaq_data.timestamps(1:2));
ind2(find((diff(ind2)./Nidaq_rate) < .01)) = [];
% remove any offsets before stim onsets - if pulse starts high
if ind2(1) < ind(1)
    ind2(1) = [];
end
frame_stim_offsets = nidaq_data.timestamps(ind2);
diff_ind2 = diff([-1*min_time_between ;frame_stim_offsets]);
ind_error2 = find(diff_ind2 < min_time_between);
frame_stim_offsets(ind_error2) = [];

if numel(frame_stim_onsets) ~= length(Condition)
    pp = warndlg(['It looks like the script has identified ' num2str(numel(frame_stim_onsets)) ...
        ' number of stimuli than is logged in the BHV file ' num2str(length(Condition)) ...
        ' Double check if your nidaq run ended before the visual stimulus ended.']);
end


%% Loop for sorting Visual Stimulus

% The struct contains all the same metadata fields as the visual
% stimulation struct:
names = fieldnames(stim.frame)';
frame_2p_metadata = cell2struct(cell(size(names)), names, 2);

% preallocate
prospect = zeros(1,n2pFrames);
for n = names
    [frame_2p_metadata(1).(n{:})] = deal(NaN(size(prospect)));
end

for f2p = 1:n2pFrames
    % Find corresponding stim frame for the current 2p frame:
	% This is defined as the last stim frame that was displayed before the
	% current 2p frame. We want to conserve causality, so it wouldn't make
	% sense to get the frame that is closest in time, because it might be
	% in the future, after the 2p frame.
	
	% If the Nidaq recording stops before the 2p recording, there are too
	% many frames in the movie. They are skipped here. The user is warned
	% about this above.
	if numel(frame_2p_onsets) < f2p
		continue
    end
    
    time_2p = frame_2p_onsets(f2p);
	frame_stim = find(frame_stim_onsets < time_2p, 1, 'last'); % Find onset timestamp prior to 2p frame
    frame_stim_offset_post = find(frame_stim_offsets > time_2p, 1, 'first');
    frame_stim_offset_last = find(frame_stim_offsets < time_2p, 1, 'last'); % Find offset prior to given 2p frame
    
    % This deals with off periods between stimuli because the last stimulus
    % onset had an offset that followed. Only want to record data when in
    % between a stimulus onset and a stimulus offset
    if frame_stim_offset_last >= frame_stim 
        continue
    end

    
    % Skip 2p frames that were acquired before stimulation started:
	if isempty(frame_stim)
		continue
    end
    
    % Skip 2p frames that were acquired after stimulation ended:
	if time_2p > frame_stim_offsets(end)
		continue
    end
    
    % Pack the stim data into the 2p metadata struct:
	for n = names
		frame_2p_metadata.(n{:})(:, f2p) = stim.frame.(n{:})(frame_stim,:);
    end	
    
    multiWaitbar('Processing frame metadata', f2p/n2pFrames);

end

% Formatting contrast so it will be compatible with makePeriRunTiffStack
frame_2p_metadata.contrast(find(isnan(frame_2p_metadata.contrast))) = 0; % Artificially setting times where stim isn't on to a value of zero
end

%% Sort Eye, Licking, Running, etc

frame_2p_metadata(1).meta.framerate_2p = 1/median(diff(frame_2p_onsets));
% save conditions and script name
frame_2p_metadata(1).meta.condition_number = unique(PPPack.hf.remove_nans(frame_2p_metadata.condition));
for i = 1:length(frame_2p_metadata(1).meta.condition_number)
    frame_2p_metadata(1).meta.condition_name{i} = stim.TimingFileByCond{...
        frame_2p_metadata(1).meta.condition_number(i)};
end




dFrame = median(diff(frame_2p_onsets));
%  want to know about specific number of licks, ensure, or quinine triggers per frame
for f2p = 1:n2pFrames; %- added RR
    if numel(frame_2p_onsets) < f2p
        frame_2p_metadata.lickingtrigsperframe(f2p) = nan;
        frame_2p_metadata.ensuretrigsperframe(f2p) = nan;
        frame_2p_metadata.quininetrigsperframe(f2p) = nan;
		continue
    end
    frame_onset = frame_2p_onsets(f2p);
    
    clear ind
    ind = find(frame_licking_onsets> (frame_onset-dFrame) & frame_licking_onsets <= (frame_onset) );
    frame_2p_metadata.lickingtrigsperframe(f2p) = length(ind);
    
    clear ind
    ind = find(frame_ensure_onsets> (frame_onset-dFrame) & frame_ensure_onsets <= (frame_onset) );
    frame_2p_metadata.ensuretrigsperframe(f2p) = length(ind);
    
    clear ind
    ind = find(frame_quinine_onsets> (frame_onset-dFrame) & frame_quinine_onsets <= (frame_onset) );
    frame_2p_metadata.quininetrigsperframe(f2p) = length(ind); 
end
multiWaitbar('Processing frame metadata', 'Close');
multiWaitbar('CloseAll');
end

function [output,ind] = get_onsets(nidaq_data,ch_to_use,thresh_to_use)
    % get the onsets for each channel type - rising edges
    if (nargin < 3)
        thresh_to_use = range(nidaq_data.data(:,ch_to_use))/2;
    end
    ind = find(diff(nidaq_data.data(:,ch_to_use)>thresh_to_use) == 1);
    Nidaq_rate = 1./diff(nidaq_data.timestamps(1:2));
    ind(find((diff(ind)./Nidaq_rate) < .01)) = [];
    output = nidaq_data.timestamps(ind);
end

% parse orientation of retinototpy script
function stim = parse_orientation_tuning_movie(stim,Condition,unique_Condition,curr_Timing_Script,count)
    ori_possibilities = [0 45 90 135 180 225 270 315 360]; % possible orientations
    for curr_num = ori_possibilities
        % 2 possible string formats
        strposs1 = sprintf('(%sd',num2str(curr_num));
        strposs2 = sprintf('Ori_%s',num2str(curr_num));
        if ~isempty(findstr(curr_Timing_Script{1},strposs1)) || ~isempty(findstr(curr_Timing_Script{1},strposs2))
            ind_curr = find(Condition == unique_Condition(count));
            stim.frame.orientation(ind_curr) = curr_num;
            stim.frame.sf(ind_curr)          = 0.04;
            stim.frame.tf(ind_curr)          = 2;
            stim.frame.location(ind_curr)    = -1;
            stim.frame.contrast(ind_curr)    = -0.8;
            break
        end
    end
end

function stim = parse_retinotopy_movies(stim,Condition,unique_Condition,curr_Timing_Script,count)
    location_possibilities = [1:6]; % possible locations
    for curr_num = location_possibilities
        strposs1 = sprintf('Pos%s',num2str(curr_num));
        strposs2 = sprintf('Pos_%s',num2str(curr_num));
         if ~isempty(findstr(curr_Timing_Script{1},strposs1)) || ~isempty(findstr(curr_Timing_Script{1},strposs2));
            ind_curr = find(Condition == unique_Condition(count));
            stim.frame.orientation(ind_curr) = 45;
            stim.frame.sf(ind_curr)          = 0.04;
            stim.frame.tf(ind_curr)          = 2;
            stim.frame.location(ind_curr)    = curr_num;
            stim.frame.contrast(ind_curr)    = -0.8;
         end
    end
end

