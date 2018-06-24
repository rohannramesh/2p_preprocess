function obj = sbxTrialInfo(obj)
% Create one variable that has all of the information on a trial by trial
% basis - each row corresponds to a single trial and each column to a
% different variable:
% Stim and Behavior Matrix = trial # x stim/behavioral variable
% ie if 100 stimuli and care about ori sf tf etc, then 100 x 14;
% Order of columns:
%     1. Stimulus Onsets (iOnsets) - frame number of these onsets
%     2. Location
%     3. Orientation
%     4. SF
%     5. TF
%     6. Ensure
%     7. Condition
%     8. Trial Error - behavior of the animal (hit (0), miss (1), correct reject (2 4), false alarm (3 5))
%     9. Block number
%     10.Contrast
%     11. Quinine
%     12. Running - speed
%     13. Pupil Diameter
%     14. Pupil Position

% iterate through each run
nRuns = length(obj.Dirs.runs);
for curr_run = 1:nRuns
    % this is the directory information for a specific run
    dirs = obj.Dirs.runs{curr_run};
    % load stim
    if ~isempty(dirs.ptb)
        load(dirs.ptb,'-mat');
    elseif ~isempty(dirs.ml) % load in MonkeyLogic file and save
        stim = PPPack.hf.bhv_read(dirs.ml);
    end
    % load .f2p file
    load(strrep(dirs.sbx,'.sbx','.f2p'),'-mat');

    % to determine if its OTB (on the ball) training
    if ~isempty(dirs.ml)
        OTB_tag = 1;
    else
        OTB_tag = 0;
    end

    frame_2p_metadata.meta.OTB_tag = OTB_tag;
    if frame_2p_metadata.meta.OTB_tag == 1
        stim.meta.frame_param_struct_settings.faststim_flag = 0;
    end
    if isfield(frame_2p_metadata,'faststim_CS') % if using special stimulus type - most users can ignore
        stim.meta.frame_param_struct_settings.faststim_flag = 1;
    end

    % now lets add running and eye if they exist
    if ~isempty(dirs.quad)
        TMP = load(dirs.quad);
        [mouse,date,run] = PPPack.hf.get_mouse_day_run_info_from_dirs(dirs);
        curr_run_speed = PPPack.hf.sbxSpeed(mouse,date,run);
        stim.RunningSpeed = curr_run_speed(1:length(frame_2p_metadata.contrast));
    end
    if ~isempty(dirs.eye_processed)
        TMP = load(dirs.eye_processed,'-mat');
        [mouse,date,run] = PPPack.hf.get_mouse_day_run_info_from_dirs(dirs);
        stim.PupilDiameter = TMP.Pupil.diameter_interp;
        stim.PupilPosition = TMP.Pupil.position_interp;
    end

    [Stim_Master,Column_Order] = PPPack.hf.getTrialInfo(frame_2p_metadata,stim,obj);


    % save
    save(strrep(dirs.sbx,'.sbx','.TrialVar'),'Stim_Master','Column_Order'); 
end

% update log file
obj.update_log_file('trial_info');
