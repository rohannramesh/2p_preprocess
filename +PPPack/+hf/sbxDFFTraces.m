function signals = sbxDFFTraces(obj,signals,frame_2p_metadata)

fr_number = length(signals.timecourse(1).raw);
nROIs = length(signals.timecourse);

%% Calculate dFF from raw timecourses
if ~isempty(frame_2p_metadata)

    fps = frame_2p_metadata.meta.framerate_2p;

    %% To Calculate f0
    if isfield(frame_2p_metadata,'contrast')
        % Identify Stimulus Onsets to get f0 window
        logic    = (frame_2p_metadata.contrast ~= 0);
        difflogic = [0 diff(logic)];
        iOnsets  = find(difflogic == 1);
        iOffsets = find(difflogic == -1);
        iOffsets(iOffsets<iOnsets(1)) = [];

        % To make iOffsets identify last frame of visual stimulus as opposed to
        % frame after visual stimulus
        iOffsets = iOffsets - 1;

        % If there's no blank period at the end of the stim period, then we have to
        % add the last offset manually:
        % (There are NaNs in the contrast field for 2p frames that were acquired
        % before or after the stimulus presentation happened. This is used to find
        % the last 2p frame that was acquired during stimulus presentation.)
        if frame_2p_metadata.contrast(end) ~= 0
            iOffsets(end+1) = ...
                find(isnan(frame_2p_metadata.contrast)==0, 1, 'last');
        end

        % How many frames before to calculate f0
        % size(signals(1).timecourse.raw,2)
        frames_before = round(fps .* obj.PreProcessingParameters.s_before);
        f0_iOnsets = iOnsets - frames_before; % identify first frame of fpre
        multiWaitbar('Processing each ROI f0', 0);
        for count = 1:nROIs % Do for each ROI
            traces.f0{count} = zeros(size(signals.timecourse(count).subtracted,1),size(signals.timecourse(count).subtracted,2));
            for iOn = 1:length(f0_iOnsets);
                % For each pre window calculate the average fluorescence 
                traces.f0{count}(1,f0_iOnsets(iOn)) = mean(signals.timecourse(count).subtracted(f0_iOnsets(iOn):f0_iOnsets(iOn)+frames_before-1));
            end
            ind_non_pre = find(traces.f0{count} == 0);

            for r = 1:(length(f0_iOnsets)-1)
                traces.f0{count}(1,f0_iOnsets(r):f0_iOnsets(r+1)-1) = traces.f0{count}(1,f0_iOnsets(r));
            end

            traces.f0{count}(1,1:f0_iOnsets(1)) = traces.f0{count}(1,f0_iOnsets(1));
            traces.f0{count}(1,f0_iOnsets(end):end) = traces.f0{count}(1,f0_iOnsets(end));
            signals.timecourse(count).f0 = traces.f0{count};
            multiWaitbar('Processing each ROI f0', count/nROIs);
        end
        multiWaitbar('Processing each ROI f0', 'Close');

        % Calculate actual dFF
        for i = 1:nROIs
            signals.timecourse(i).dff = (signals.timecourse(i).subtracted - signals.timecourse(i).f0)./(signals.timecourse(i).f0);
            signals.timecourse(i).dff_norm = signals.timecourse(i).dff ./ max(signals.timecourse(i).dff);
        end
    end
end
%% Now calculate dFF using sliding window method

% Calculate f0 for each timecourse using a moving window of time window
% prior to each frame
f0_vector = zeros(nROIs,fr_number);
time_window = 32; % moving window of X seconds - calculate f0 at time window prior to each frame - used to be 32
percentile = 10; % used to be 30
time_window_frame = round(time_window*frame_2p_metadata.meta.framerate_2p);

% create temporary traces variable that allows you to do the prctile
% quickly
traces_f = nan(nROIs,length(signals.timecourse(1).subtracted));
for curr_ROI = 1:nROIs
    traces_f(curr_ROI,:) = signals.timecourse(curr_ROI).subtracted;
end

% parallel processing
try a = gcp('nocreate');
catch err
    a = [];
end
if isempty(a)
    poolobj = parpool(10);
    pool_siz = poolobj.NumWorkers;;
else
    pool_siz = a.NumWorkers;
    poolobj = a;
end

% how many ROIs per core
nROIs_per_core = ceil(nROIs/pool_siz);
ROI_vec = 1:nROIs_per_core.*pool_siz;
ROI_blocks = PPPack.hf.unshuffle_array(ROI_vec,nROIs_per_core);
ROI_start_points = ROI_blocks(:,1);   
parfor curr_ROI_ind = 1:pool_siz
    ROIs_to_use = ROI_blocks(curr_ROI_ind,:);
    ROIs_to_use(ROIs_to_use > nROIs) = [];
    % pre-allocate
    f0_vector_cell{curr_ROI_ind} = nan(length(ROIs_to_use),fr_number);
    for i = 1:fr_number
        if i <= time_window_frame
            frames = traces_f(ROIs_to_use,1:time_window_frame);
            f0 = prctile(frames,percentile,2);
        else
            frames = traces_f(ROIs_to_use,i - time_window_frame:i-1);
            f0 = prctile(frames,percentile,2);
        end
        f0_vector_cell{curr_ROI_ind}(:,i) = f0;
    end
end
% reshape into correct structure
for curr_ROI_ind = 1:pool_siz
    ROIs_to_use = ROI_blocks(curr_ROI_ind,:)
    ROIs_to_use(ROIs_to_use > nROIs) = [];
    f0_vector(ROIs_to_use,:) = f0_vector_cell{curr_ROI_ind};
end


traces_dff = (traces_f-f0_vector)./ f0_vector;

%stick back into signals variable
for curr_ROI = 1:nROIs
    signals.timecourse(curr_ROI).f0_axon = f0_vector(curr_ROI,:);
    signals.timecourse(curr_ROI).dff_sliding = traces_dff(curr_ROI,:);
    signals.timecourse(curr_ROI).dff_sliding_norm =  signals.timecourse(curr_ROI).dff_sliding ./ ...
        max( signals.timecourse(curr_ROI).dff_sliding);
end    

if ~isempty(poolobj)
    delete(poolobj);
end
