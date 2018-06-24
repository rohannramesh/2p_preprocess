function [Stim_Master,Column_Order] = getTrialInfo(frame_2p_metadata,stim,obj)
%% Create Master Matrices with all of Stim and Behavioral Information & dFF information

% Stim and Behavior Matrix = trial # x stim/behavioral variable
% ie if 100 stimuli and care about ori sf tf etc, then 100 x 5;
% Order of columns:
%     1. Stimulus Onsets (iOnsets)
%     2. Location
%     3. Orientation
%     4. SF
%     5. TF
%     6. Ensure
%     7. Condition
%     8. Trial Error
%     9. Block number
%     10.Contrast
%     11. Quinine
%     12. Running
%     13. Pupil Diameter
%     14. Pupil Position

Column_Order = {'iOnsets','Location','Orientation','SF','TF','Ensure','Condition','TrialError',...
    'Block Number','Contrast','Quinine','Running','Pupil Diameter','Pupil Position'};

        
if isfield(frame_2p_metadata.meta, 'OTB_tag')
    OTB_tag = frame_2p_metadata.meta.OTB_tag;
else
    OTB_tag = 0;
end

column_number = 14;

% find each stimulus onset - this variable will be called iOnsets
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

Stim_Master = nan(length(iOnsets),column_number);

% Stimulus Onsets is column 1
Stim_Master(:,1) = iOnsets;

% Identify Location at each stimulus - using Location_Sum
if size(frame_2p_metadata.location,1) > 1 % Added because OTB training has one vector for location
    frame_2p_metadata.location_sum = sum(frame_2p_metadata.location);
else
    frame_2p_metadata.location_sum = frame_2p_metadata.location;
end
Stim_Master(:,2) = frame_2p_metadata.location_sum(iOnsets);

% Identify Orientation at each stimulus
Stim_Master(:,3) = frame_2p_metadata.orientation(iOnsets);

% Identify SF at each stimulus
Stim_Master(:,4) = frame_2p_metadata.sf(iOnsets);

% Identify TF at each stimulus
try Stim_Master(:,5) = frame_2p_metadata.temporal_frequencies_hz(iOnsets);
catch err
    display('Assuming it is OTB training')
    Stim_Master(:,5) = frame_2p_metadata.tf(iOnsets);
end

% Identify OTB variable that correspond to each iOnset
if OTB_tag == 1;
    % Identify Ensure onset for each stim if exists - column 6
    Ensure_Onsets = find(frame_2p_metadata.ensuretrigsperframe == 1);
    for count = 1:length(Ensure_Onsets);
        curr_onset = Ensure_Onsets(count);
        ind_place_val = find(iOnsets < curr_onset,1,'last');
        Stim_Master(ind_place_val,6) = curr_onset;
    end
    % Identify condition present for each stimulus - number that
    % corresponds to trial type
    Stim_Master(:,7) = frame_2p_metadata.condition(iOnsets);
    % Identify Trial Error - behavior of the animal
    if isfield(frame_2p_metadata,'trialerror')
        Stim_Master(:,8) = frame_2p_metadata.trialerror(iOnsets);
    else
        Stim_Master(:,8) = stim.TrialError(1:length(iOnsets));
    end
    
    % Identify Block Number
    if isfield(frame_2p_metadata,'block')
        Stim_Master(:,9) = frame_2p_metadata.block(iOnsets);
    else
        Stim_Master(:,9) = stim.BlockNumber(1:length(iOnsets));
    end
    
    % Contrast
    Stim_Master(:,10) = frame_2p_metadata.contrast(iOnsets);
    
    % quinine
    if isfield(frame_2p_metadata,'quininetrigsperframe')
            % Identify Ensure onset for each stim if exists - column 6
            Quinine_Onsets = find(frame_2p_metadata.quininetrigsperframe == 1);
            for count = 1:length(Quinine_Onsets);
                curr_onset = Quinine_Onsets(count);
                ind_place_val = find(iOnsets < curr_onset,1,'last');
                Stim_Master(ind_place_val,11) = curr_onset;
            end
    else
        Stim_Master(:,11) = nan(size(iOnsets));
    end
    
    % running 
    fps = frame_2p_metadata.meta.framerate_2p;
    Onesec_pre_onsets = iOnsets - round(fps*obj.PreProcessingParameters.s_before);
    Onesec_pre_onsets(Onesec_pre_onsets < 0) = 1; % don't let there be neg values
    if isfield(stim,'RunningSpeed')
        Running_responses = nan(1,length(iOnsets));
        for ii = 1:length(Onesec_pre_onsets)
            Running_responses(ii) = nanmean(stim.RunningSpeed(Onesec_pre_onsets(ii):iOnsets(ii)-1));
        end
        Stim_Master(:,12) = Running_responses;
    else
        Stim_Master(:,12) = nan(size(iOnsets));
    end
    
    % Pupil Diameter
    if isfield(stim,'PupilDiameter')
        PupilDiameter = nan(1,length(iOnsets));
        for ii = 1:length(Onesec_pre_onsets)
            PupilDiameter(ii) = nanmean(stim.PupilDiameter(Onesec_pre_onsets(ii):iOnsets(ii)-1));
        end
        Stim_Master(:,13) = PupilDiameter;
    else
        Stim_Master(:,13) = nan(size(iOnsets));
    end
    
    % Pupil Position
    if isfield(stim,'PupilPosition')
        PupilPosition = nan(1,length(iOnsets));
        for ii = 1:length(Onesec_pre_onsets)
            PupilPosition(ii) = nanmean(stim.PupilPosition(Onesec_pre_onsets(ii):iOnsets(ii)-1));
        end
        Stim_Master(:,14) = PupilPosition;
    else
        Stim_Master(:,14) = nan(size(iOnsets));
    end            
else
    % everything except contrast is behavioral variables
    Stim_Master(:,6) = nan(size(iOnsets));
    Stim_Master(:,7) = nan(size(iOnsets));
    Stim_Master(:,8) = nan(size(iOnsets));
    Stim_Master(:,9) = nan(size(iOnsets));
    % contrast
    Stim_Master(:,10) = frame_2p_metadata.contrast(iOnsets);
    Stim_Master(:,11) = nan(size(iOnsets));
    % pupil and running
    Stim_Master(:,12) = nan(size(iOnsets));
    Stim_Master(:,13) = nan(size(iOnsets));
    Stim_Master(:,14) = nan(size(iOnsets));
end




