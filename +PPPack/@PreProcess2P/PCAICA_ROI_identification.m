function obj = PCAICA_ROI_identification(obj)
% Run Mukamel and Schnitzer's PCA/ICA ROI identification algorithm
% This algorithm will identify putative ROIs after de-noising the data via
% PCA and ROI identification using ICA (implemented by Rohan Ramesh, see
% Mukamel and Schnitzer, 2009) 
% Output of running this will save a *.ica file in the directory of the
% last run included. This file will be formatted to either be run through
% the Matlab GUI, preprocessed for the Web Clicking GUI, or run through a
% pre-trained Convolutional Neural Network for ROI selection

% load in mov from all runs
mov_ds = []; % this will be the across runs mov variable
% iterate through each run
nRuns = length(obj.Dirs.runs);
for curr_run = 1:nRuns
    % this is the directory information for a specific run
    curr_dir = obj.Dirs.runs{curr_run};
    
    % for now also cd into the proper directory
    cd(curr_dir.path)

    % load in mov - crop edges, spatially downsample by a factor of 2, and
    % temporally downsample by a factor of 5
    info = PPPack.hf.sbxInfo(curr_dir.sbx);
    tic
    mov =  PPPack.hf.sbxLoadDSMov(obj,curr_dir);
    toc
    sz_mov = size(mov);
    mov_ds = cat(3,mov_ds,mov);
end


% now get cell masks using mukamel and schnitzer
tempmovpath = fullfile(curr_dir.path, 'temp.tif');
if isfield(obj.PreProcessingParameters.ROI_type,'ROI_type')
    if strcmp(obj.PreProcessingParameters.ROI_type,'Cell')
        PPPack.PCAICA.run_PCAICA(tempmovpath,curr_dir.path,mov_ds);
    elseif strcmp(obj.PreProcessingParameters.ROI_type,'Axon')
        PPPack.PCAICA.run_PCAICA_Axon(tempmovpath,curr_dir.path,mov_ds);
    end
else % default assumption is Cell
    PPPack.PCAICA.run_PCAICA(tempmovpath,curr_dir.path,mov_ds);
end
% Save icaguidata for Matlab GUI
global premasks
save([curr_dir.sbx_name '.ica'], 'premasks','-v7.3');

% preprocess for javascript
[mouse,date,run] = PPPack.hf.get_mouse_day_run_info_from_dirs(curr_dir);
PPPack.hf.processForJavascript(mouse,date,run,[],[],obj.PreProcessingParameters.server);

% update log file
obj.update_log_file('PCAICA_ROI_identification');