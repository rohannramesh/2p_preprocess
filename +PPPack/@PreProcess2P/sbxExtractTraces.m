function obj = sbxExtractTraces(obj)
% lets load in the premasks file whether from PCA/ICA or from CNMF.
% Assuming you have either cell clicked using the custom Matlab GUI, the
% CNN, or the javascript GUI and the .ica or the .nmf variables have been
% made - to make these see PCAICA_ROI_identification or CNMF_ROI_identification
% Then extract fluorescense timecourses from the sbx file, do neuropil
% corrections (either weighted or non-weighted), and calculate the dFF
% timecourse (both sliding window and re-zeroed prior to visual stimulus
% onset - if exists) and save

% first lets check which type: PCA/ICA vs NMF and load in premasks
clearvars -global
if strcmp(obj.PreProcessingParameters.ROI_algorithm,'PCA/ICA')
    display('PCA/ICA extraction')
    % load premasks
    premasks_filename = strrep(obj.Dirs.runs{end}.sbx,'.sbx','.icamasks');
    load(premasks_filename,'-mat');
elseif strcmp(obj.PreProcessingParameters.ROI_algorithm,'NMF')
    display('NMF extraction')
    namefile = dir([obj.Dirs.date_mouse '*.nmfmasks']);
    premasks_filename = [obj.Dirs.date_mouse namefile.name];
    load(premasks_filename,'-mat');
else
    error('Incorrect ROI algorithm')
end


%% Extract

% iterate through each run
nRuns = length(obj.Dirs.runs);
for curr_run = 1:nRuns
    % this is the directory information for a specific run
    dirs = obj.Dirs.runs{curr_run};
    
    % load f2p variable
    try load(strrep(dirs.sbx,'.sbx','.f2p'),'-mat')
    catch err
        frame_2p_metadata = [];
        info = PPPack.hf.sbxInfo(dirs.sbx,1);
        if info.scanmode == 0
            frame_2p_metadata.meta.framerate_2p = 30.96;
            warning('You didnt save nidaq data so checking the info variable and hardcoding 31')
        elseif info.scanmode == 1
            frame_2p_metadata.meta.framerate_2p = 15.48;
            warning('You didnt save nidaq data so checking the info variable and hardcoding 15')
        end
    end
    
    % extract traces
    [sig,signals] = PPPack.hf.sbxGetRawTraces(obj,dirs,premasks)

    % calculate dFF
    signals = PPPack.hf.sbxDFFTraces(obj,signals,frame_2p_metadata)
    
    % try doing deconvolution
    try signals = PPPack.hf.deconvolveRemovingNANRows(signals);
    catch err
        warning('make sure initialized cvx etc correctly for NMF')
    end
    
    % make image of FOV plus ROIs
    PPPack.hf.plot_overlay_image_masks(dirs,obj,signals);
    
    % save
    if strcmp(obj.PreProcessingParameters.ROI_algorithm,'PCA/ICA')
        save(strrep(dirs.sbx,'.sbx','.signalsica'),'signals','-v7.3');
    elseif strcmp(obj.PreProcessingParameters.ROI_algorithm,'NMF')
        save(strrep(dirs.sbx,'.sbx','.signalsnmf'),'signals','-v7.3');
    end
end

% update log file
obj.update_log_file('trace_extraction');