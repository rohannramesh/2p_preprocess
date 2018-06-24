function obj = sbxCreateMetadata_dFFimages(obj)
% if running the Monkeylogic behavioral set up (commonly used for training
% in the Andermann lab) or using a Psychtoolbox visual stimulus set up this
% puts the visual stimulus and behavior information into the desired 2P
% framerate. 
% This script will have to be adjusted across labs to deal with different
% behavioral paradigms

% iterate through each run
nRuns = length(obj.Dirs.runs);
for curr_run = 1:nRuns
    % this is the directory information for a specific run
    dirs = obj.Dirs.runs{curr_run};
    
    % only run if there is an ephys (nidaq) file type and either
    % psychtoolbox (ptb) or monkeylogic (ml) file
    if ~isempty(dirs.nidaq) && (~isempty(dirs.ptb) || ~isempty(dirs.ml)) % if there are pulses and a stim file
        % get the info file for that run
        info = PPPack.hf.sbxInfo(dirs.sbx);
    end