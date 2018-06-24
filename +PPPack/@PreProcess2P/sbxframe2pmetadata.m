function obj = sbxframe2pmetadata(obj)
% Make the metadata variable that contains information about the visual
% information the animal is seeing on any given 2p frame as well as what
% the animal is doing behaviorally on any given frame (licking) - running
% and eye tracking are handled separately due to file type
% currently integrates animal behavior using the Monkeylogic system (Asaad
% and Eskandar)

% iterate through each run
nRuns = length(obj.Dirs.runs);
for curr_run = 1:nRuns
    % this is the directory information for a specific run
    dirs = obj.Dirs.runs{curr_run};

    % if frame_2p_metadata already exists then don't load it in again
    if exist([dirs.path '\' dirs.sbx_name '.f2p'])
        display('Already found frame_2p_metadata in folder')
        continue
    end

    % lets get the size of the movie
    info = PPPack.hf.sbxInfo(dirs.sbx);
    szmov = info.max_idx+1;

    % Deal with difference in formatting - 
    % nidaq files if you are using a separately designed computer to align 2p
    % pulses and stimulus/ behavioral information
    % ephys files if you are using the system implemented withing the Scanbox
    % software
    if ~isempty(dirs.nidaq) && strcmp(dirs.nidaq(end-4:end),'nidaq')
        nidaq = load(dirs.nidaq,'-mat');
    elseif ~isempty(dirs.nidaq) && strcmp(dirs.nidaq(end-4:end),'ephys')
        nidaq = PPPack.hf.process_ephys_files(dirs.nidaq);
    else
        frame_2p_metadata = [];
        return
    end
    % just for buffering purposes in the get2pFrameMetadata script will buffer
    % with zeros
    sz = size(nidaq.data);
    if sz(2) < 8
        stuff_add = 8-sz(2);
        nidaq.data(:,8-stuff_add+1:8) = 0;
    end
    nidaq.timestamps = nidaq.timeStamps;
    nidaq = rmfield(nidaq,'timeStamps');
    % these are your channels of interest
    RigName = obj.PreProcessingParameters.rig
    ch_values = PPPack.hf.get_ch_values_nidaq(RigName);
    % try loading in quad data and if not empty pass to following functions in
    % nidaq
    if ~isempty(dirs.quad)
        TMP = load(dirs.quad);
        nidaq.running = TMP.quad_data;
    end
    % load in the stimulus information
    if ~isempty(dirs.ml) % load in MonkeyLogic file and save
        stim = PPPack.hf.bhv_read(dirs.ml);
        % now need to determine what type of run it is - Conditions, Ori, or Retinotopy
        if length(strfind(stim.ConditionsFile,'Gabor_Orientation')) > 0
            stim.RunType = 'Ori'
        elseif length(strfind(stim.ConditionsFile,'Retinotopy')) > 0
            stim.RunType = 'Retinotopy'
        else
            stim.RunType = 'Training'
        end
        % get frame_2p_metadata for Monkeylogic
        frame_2p_metadata = PPPack.hf.get2pFrameMetadata_ml(stim, nidaq, szmov,ch_values);
    elseif ~isempty(dirs.ptb) % if visual stimulus is a pyschtoolbox file
        load(dirs.ptb,'-mat');
        % get frame_2p_metadata for PTB
        error('Script needs to be optimized for PTB before running in this package')
    %     frame_2p_metadata = PPPack.hf.get2pFrameMetadata_ptb(stim, nidaq, szmov,[],ch_values);
    else
        frame_2p_metadata = [];
        return
    end


    % if people have screwed up the contrast and sf - this will check and fix
    % these metrics for later sorting purposes
    unique_sf = unique(frame_2p_metadata.sf(~isnan(frame_2p_metadata.sf)));
    try unique_tf = unique(frame_2p_metadata.temporal_frequencies_hz(~isnan(frame_2p_metadata.temporal_frequencies_hz)));
    catch err
        unique_tf = unique(frame_2p_metadata.tf(~isnan(frame_2p_metadata.tf)));
    end
    % if one of the sf = 0 something has gone wrong in your stim script and
    % correcting now
    if sum(unique_sf == 0) > 0
        ind_error = find(frame_2p_metadata.orientation == 360);
        frame_2p_metadata.sf(ind_error) = unique_sf(2);
        try frame_2p_metadata.temporal_frequencies_hz(ind_error) = unique_tf(2);
        catch err
            frame_2p_metadata.tf(ind_error) = unique_tf(2);
        end
        warning('Correcting your blank trials bc your stim script is messed up. Dont set the sf and tf to 0. Assuming only 1 sf and tf')
    end

    save([dirs.path '\' dirs.sbx_name '.f2p'],'frame_2p_metadata')
end

% update log file
obj.update_log_file('frame_2p_metadata');
