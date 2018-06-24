function frame_2p_metadata = get2pFrameMetadata_ptb(stim, nidaq_data, n2pFrames, eye_data,ch)
% get2pFrameMetadata takes stimulation metadata and NiDAQ-recorded timing
% data for a 2p movie and returns a struct with metadata for each frame in
% the movie.
%
% INPUTS:
%	stim		Visual stimulation metadata struct as saved by Andermann
%               lab stimulation script.
%	nidaq_data	Timing data recorded via NiDAQ during 2p acquisition.
%               For channel number see get_ch_values
%
%	n2pFrames	The recorded number of 2p frames. The number of frame
%				pulses that appear in the NiDAQ data are a few more than 
%				the recorded number of frames. Testing showed that these
%				extra pulses occur at the end of 2p acquisition such that
%				the first recorded pulse corresponds to the first recorded
%				frame.

%	
% OUTPUTS:
%	frame_2p_metadata	The structure that contains all the metadata for
%						the 2p frames.



%% Sort Nidaq data
% The Nidaq data -- unbelievably -- is not saved in chronological order by
% the recording device. This is not an issue with regards to accuracy
% because timestamps are saved (correctly) for each data point. However, 
% We rely on chronological order below, so we have to sort the data here.
% This does not change the timing of each data point, it just changes the 
% order in which the data are saved in the arrays.
try [nidaq_data.timestamps, time_idx] = sort(nidaq_data.timestamps);

    for i = 1:size(nidaq_data.data, 2)
        nidaq_data.data(:, i) = nidaq_data.data(time_idx, i);
    end
catch err
end

%% Extract running eyetracking and bodytracking events
thresh_running = range(nidaq_data.data(:,ch.running))/2;
ind = find(diff(nidaq_data.data(:,ch.running)>thresh_running) == 1);
try Nidaq_rate = 1./diff(nidaq_data.timestamps(1:2));
catch err
    Nidaq_rate = nidaq_data.recording_frequency;
end
ind(find((diff(ind)./Nidaq_rate) < .01)) = [];
frame_running_onsets = nidaq_data.timestamps(ind);

thresh_eyetrack = range(nidaq_data.data(:,ch.eyetrack))/2;
ind = find(diff(nidaq_data.data(:,ch.eyetrack)>thresh_eyetrack) == 1);
% Nidaq_rate = 1./diff(nidaq_data.timestamps(1:2));
ind(find((diff(ind)./Nidaq_rate) < .01)) = [];
frame_eyetrack_onsets = nidaq_data.timestamps(ind);
frame_eyetrack_pulses = zeros(length(nidaq_data.data(:,ch.eyetrack)),1);
frame_eyetrack_pulses(ind) = 1;

clear ind
thresh_licking = range(nidaq_data.data(:,ch.licking))/2;
ind = find(diff(nidaq_data.data(:,ch.licking)>thresh_licking) == 1);
% Nidaq_rate = 1./diff(nidaq_data.timestamps(1:2));
ind(find((diff(ind)./Nidaq_rate) < .01)) = [];
frame_licking_onsets = nidaq_data.timestamps(ind);

clear ind
thresh_ensure = range(nidaq_data.data(:,ch.ensure))/2;
ind = find(diff(nidaq_data.data(:,ch.ensure)>thresh_ensure) == 1);
% Nidaq_rate = 1./diff(nidaq_data.timestamps(1:2));
ind(find((diff(ind)./Nidaq_rate) < .01)) = [];
frame_ensure_onsets = nidaq_data.timestamps(ind);

% % ONLY USE THIS WHEN REFERENCING BODYTRACKING
% ind_errors = find(diff(frame_eyetrack_onsets)>0.1); % find instances where gap between onsets are large
% ind_errors_max_start = ind_errors((max(find(ind_errors < 1000))));
% ind_errors_min_end = ind_errors((min(find(ind_errors > length(frame_eyetrack_onsets)-1000))));
% 
% 
% frame_eyetrack_onsets(1:ind_errors_max_start) = [];
% ind(1:ind_errors_max_start) = [];
% 
% frame_eyetrack_onsets(ind_errors_min_end-1:end) = [];
% ind(ind_errors_min_end-1:end) = [];
% frame_eyetrack_pulses_corrected = zeros(length(nidaq_data.data(:,ch.eyetrack)),1);
% frame_eyetrack_pulses_corrected(ind) = 1;
% if ~isempty(eye_data)
%     figure;
%     ax(1) = subplot(2,1,1); plot(frame_eyetrack_pulses)
%     ax(2) = subplot(2,1,2); plot(frame_eyetrack_pulses_corrected)
%     linkaxes(ax, 'x')
% end

thresh_bodytrack = range(nidaq_data.data(:,ch.bodytrack))/2;
ind = find(diff(nidaq_data.data(:,ch.bodytrack)>thresh_bodytrack) == 1);
% Nidaq_rate = 1./diff(nidaq_data.timestamps(1:2));
ind(find((diff(ind)./Nidaq_rate) < .01)) = [];
frame_bodytrack_onsets = nidaq_data.timestamps(ind);

%% Extract times of stimulation frames (monitor framerate) from rec_data:

% Use simple thresholding to get frame onset times from the analog frame
% signal recorded by the NiDAQ system:
clear ind
thresh_2p = range(nidaq_data.data(:,ch.twoP))/2;
ind = find(diff(nidaq_data.data(:,ch.twoP)>thresh_2p) == 1);
% Nidaq_rate = 1./diff(nidaq_data.timestamps(1:2));
ind(find((diff(ind)./Nidaq_rate) < .01)) = [];
frame_2p_onsets = nidaq_data.timestamps(ind);

% if the number of f2p onsets is half the number of total movie frames
% this is bc the new 2p rig sends out only one pulse
if length(frame_2p_onsets)*2 == n2pFrames | (length(frame_2p_onsets)-1)*2 == n2pFrames % need or bc sometimes ends in up state
    frame_2p_onsets = interp1(1:2:length(frame_2p_onsets)*2,frame_2p_onsets,1:length(frame_2p_onsets)*2,'linear','extrap');
    warning('Assuming on new 2p rig. Bc only send out one pulse for forward and back scan, artificially adding pulse for flyback.')
end
    

if numel(frame_2p_onsets) < n2pFrames
	w = warndlg(['Number of 2p onsets ' num2str(length(frame_2p_onsets)) ...
        ' number of frames in movie ' num2str(n2pFrames)]);
% 	waitfor(w);
end

if ~isempty(stim)
% Visstim pulses are +1V for all frames and +5V for non-blank frames
clear ind
thresh_stim = 0.5; % Get all frames by thresholding at 0.5 V
ind = find(diff(nidaq_data.data(:,ch.visstim)>thresh_stim) == 1);
% Nidaq_rate = 1./diff(nidaq_data.timestamps(1:2));
ind(find((diff(ind)./Nidaq_rate) < .01)) = [];
frame_stim_onsets = nidaq_data.timestamps(ind);

if numel(frame_stim_onsets) == 1 + numel(stim.frame.contrast)
	warning(['For a while, a bug in the stimulation code sent out an ',...
				'additional frame pulse before the start of stimulation. ',...
				'This seems to be the case here and has been corrected automatically.']);
elseif numel(frame_stim_onsets) ~= numel(stim.frame.contrast)
	s = sprintf([...
		'The number of frame pulses detected in the Nidaq recording (%d) ',...
		'does not match the number of frames saved in the stimulation ',...
		'metadata structure (%d). This needs to be checked unless it is ',...
		'expected, e.g. due to premature Nidaq recording ending.'], ...
		numel(frame_stim_onsets),...
		numel(stim.frame.contrast));
	w = warndlg(s);
% 	waitfor(w);
end

% frame_stim_onsets = frame_stim_onsets(...
% 	diff([frame_stim_onsets; 0]) < ...
% 	0.5 .* median(diff(stim.frame.high_prec_sys_ts)));

frame_stim_interval = median(diff(frame_stim_onsets));

%% Sort all metadata into a struct:
% The struct contains all the same metadata fields as the visual
% stimulation struct:
names = fieldnames(stim.frame)';
%frame_2p_metadata = struct(names{:}, {});
frame_2p_metadata = cell2struct(cell(size(names)), names, 2);



% Clean up "location" paramerter. This was unfortunately originally
% specified with the frame number in the first rather than second
% dimension, so we change that here to adapt it to the convention of all
% other parameters:
if isfield(stim.frame, 'location') && size(stim.frame.location, 2) == 2
	stim.frame.location = stim.frame.location';
end

% Pre-allocate. Elements without metadata will be NaN:
for n = names
		[frame_2p_metadata(1).(n{:})] = deal(NaN(size(stim.frame.(n{:}))));
end

% For each 2p frame, search and store corresponding visual stimulation
% metadata:
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
	frame_stim = find(frame_stim_onsets < time_2p, 1, 'last');
	
	% Find stimulation frame index: Get stim-PC frame timestamp that is
	% closest to corresponding 2p frame (using the NiDAQ recording as a
	% common time base)
% 	frame_stim = find(frame_stim_onsets - time_2p) == ...
% 		min(abs(frame_stim_onsets - time_2p)), 1);
	
	% Skip 2p frames that were acquired before stimulation started:
	if isempty(frame_stim)
		continue
	end
	
	% Skip 2p frames that were acquired after stimulation ended:
	if time_2p > frame_stim_onsets(end) + frame_stim_interval
		continue
	end
	

	% Pack the stim data into the 2p metadata struct:
	for n = names
		frame_2p_metadata.(n{:})(:, f2p) = stim.frame.(n{:})(:, frame_stim);
	end	
	
	multiWaitbar('Processing frame metadata', f2p/n2pFrames);
end
end

frame_2p_metadata(1).meta.framerate_2p = 1/median(diff(frame_2p_onsets));


if ~isempty(eye_data)
    temporal_downsampling = eye_data.NdownT;
    ALL_mat = {eye_data.Horiz_eyepos00 eye_data.Vert_eyepos00...
        eye_data.PupRadius00};
    length_ALL = length(ALL_mat);
    for n = 1:length(ALL_mat)
        new_mat_index = n+length_ALL;
        ALL_mat{new_mat_index} = medfilt_RR(ALL_mat{n},1,15,temporal_downsampling);
    end

    eye_data.Horiz_eyepos00_medfilt = ALL_mat{4};
    eye_data.Vert_eyepos00_medfilt = ALL_mat{5};
    eye_data.PupRadius00_medfilt = ALL_mat{6};

    % Interpolate the downsampled eyetracking data set to every other frame not
    % analyzed
    clear ALL_mat
    ALL_mat = {eye_data.Horiz_eyepos00_medfilt eye_data.Vert_eyepos00_medfilt...
        eye_data.PupRadius00_medfilt double(eye_data.Mask_blinkorother00)};
    length_ALL = length(ALL_mat);
    if temporal_downsampling == 2;
        for i = 1:length(ALL_mat)
            new_mat_index = i+length_ALL;
            desired_length = length(ALL_mat{i})*2;
            x = 1:2:desired_length;
            xq = 1:desired_length;
            ALL_mat{new_mat_index} = interp1(x,ALL_mat{i},xq);
        %     figure;plot(x,ALL_mat{i},'o',xq,ALL_mat{new_mat_index});
        end

        % For logical mask_blink i.e. what problem points to be aware of, want in
        % between values to be zero
        % length_downs = length(eye_data.Mask_blinkorother00);
        % eye_data.Mask_blinkorother00_undownsampled = zeros(length_downs*2,1);
        % eye_data.Mask_blinkorother00_undownsampled(1:2:length_downs*2) =  eye_data.Mask_blinkorother00;


        eye_data.Horiz_eyepos00_undownsampled = ALL_mat{5}';
        eye_data.Vert_eyepos00_undownsampled = ALL_mat{6}';
        eye_data.PupRadius00_undownsampled = ALL_mat{7}';
        eye_data.Mask_blinkorother00_undownsampled = (ALL_mat{8} == 1)'; % Just connecting the NaNs
    else
        eye_data.Horiz_eyepos00_undownsampled = eye_data.Horiz_eyepos00_medfilt;
        eye_data.Vert_eyepos00_undownsampled = eye_data.Vert_eyepos00_medfilt;
        eye_data.PupRadius00_undownsampled = eye_data.PupRadius00_medfilt;
        eye_data.Mask_blinkorother00_undownsampled = eye_data.Mask_blinkorother00'; % Just connecting the NaNs
    end

end

%added MA
dFrame = median(diff(frame_2p_onsets));
% for f2p = 1:n2pFrames
for f2p = 1:numel(frame_2p_onsets); %- added RR
    frame_onset = frame_2p_onsets(f2p);

    ind = find(frame_running_onsets> (frame_onset-dFrame) & frame_running_onsets <= (frame_onset) );
    frame_2p_metadata.runningtrigsperframe(f2p) = length(ind);
    
    clear ind
    ind = find(frame_licking_onsets> (frame_onset-dFrame) & frame_licking_onsets <= (frame_onset) );
    frame_2p_metadata.lickingtrigsperframe(f2p) = length(ind);
    
    clear ind
    ind = find(frame_ensure_onsets> (frame_onset-dFrame) & frame_ensure_onsets <= (frame_onset) );
    frame_2p_metadata.ensuretrigsperframe(f2p) = length(ind);
    
    %nw get time of first and last movie trigger before 2p framnetrigger, for eyetrakcing
    %and bodytracking

    ind = find(frame_eyetrack_onsets> (frame_onset-dFrame) & frame_eyetrack_onsets <= (frame_onset) );
%     if ~isempty(eye_data)
%          if ind > length(eye_data.Horiz_eyepos00_undownsampled) % to get around bug
%              disp('Have more eyetrack pulses than actual frames, assuming clipped from end')
%              continue
%          end
%     end
    if ~isempty(ind)
        frame_2p_metadata.eyetrackframe_last(f2p) = frame_eyetrack_onsets(ind(end));
        frame_2p_metadata.eyetrackframe_first(f2p) = frame_eyetrack_onsets(ind(1));
        if ~isempty(eye_data)
            frame_2p_metadata.Horiz_eyepos(f2p) = mean(eye_data.Horiz_eyepos00_undownsampled(ind));
            frame_2p_metadata.Vert_eyepos(f2p) = mean(eye_data.Vert_eyepos00_undownsampled(ind));
            frame_2p_metadata.PupRadius(f2p) = mean(eye_data.PupRadius00_undownsampled(ind));
            frame_2p_metadata.Mask_blinkorother(f2p) = max(eye_data.Mask_blinkorother00_undownsampled(ind));
        end
    end
    
    ind = find(frame_bodytrack_onsets> (frame_onset-dFrame) & frame_bodytrack_onsets <= (frame_onset) );
    if ~isempty(ind)
        frame_2p_metadata.bodytrackframe_last(f2p) = frame_bodytrack_onsets(ind(end));
        frame_2p_metadata.bodytrackframe_first(f2p) = frame_bodytrack_onsets(ind(1));
    end
    
end

frame_2p_metadata.frame_bodytrack_onsets = frame_bodytrack_onsets;
frame_2p_metadata.frame_eyetrack_onsets = frame_eyetrack_onsets;

if ~isempty(eye_data)
    ALL_mat2 = {frame_2p_metadata.Horiz_eyepos frame_2p_metadata.Vert_eyepos ...
        frame_2p_metadata.PupRadius};

    if ~isempty(eye_data)
        length_ALL2 = length(ALL_mat2);
        for i = 1:length(ALL_mat2)
            new_mat_index = i+length_ALL2;
            len = length(ALL_mat2{i});
            nonzero_ind = find(ALL_mat2{i}>0);
            ind_nan = find(isnan(ALL_mat2{i}));
            xq = 1:len;
            ALL_mat2{new_mat_index} = interp1(nonzero_ind,ALL_mat2{i}(nonzero_ind),xq);
        %     ALL_mat2{new_mat_index}(ind_nan) = NaN;
            figure;plot(nonzero_ind,ALL_mat2{i}(nonzero_ind),'o',xq,ALL_mat2{new_mat_index})
            title('eyetracking data blue circles - median filtered and interpolated')
        end

        ind_nan = find(frame_2p_metadata.Mask_blinkorother == 1);
        frame_2p_metadata.Horiz_eyepos = ALL_mat2{4};
        frame_2p_metadata.Horiz_eyepos(ind_nan) = NaN; % Since Interpolated - NaN disappear - this redefines them
        frame_2p_metadata.Vert_eyepos = ALL_mat2{5};
        frame_2p_metadata.Vert_eyepos(ind_nan) = NaN;
        frame_2p_metadata.PupRadius = ALL_mat2{6};
        frame_2p_metadata.PupRadius(ind_nan) = NaN;
    end
end

multiWaitbar('Processing frame metadata', 'Close');
multiWaitbar('CloseAll');

% Making it so that the baseline for orientation is nan and not  0:
if isfield(frame_2p_metadata, 'orientation')
%     ind_stim_on = find(frame_2p_metadata.contrast == -0.8);
    ind_stim_on = find(frame_2p_metadata.contrast ~= 0);
    if isempty(ind_stim_on)
        ind_stim_on = find(frame_2p_metadata.contrast == 0.8);
    end
    tmp_ori = nan(size(frame_2p_metadata.orientation));
    tmp_ori(ind_stim_on) = frame_2p_metadata.orientation(ind_stim_on);
    frame_2p_metadata.orientation = tmp_ori;
end



