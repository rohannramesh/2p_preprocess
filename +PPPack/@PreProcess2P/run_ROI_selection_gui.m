function obj = run_ROI_selection_gui(obj)
% load in premasks variable and run either the PCAICA gui or the CNMF gui
% if have already run through the CNN uses these ROIs as suggestions for
% those to keep

% first lets check which type: PCA/ICA vs NMF
clearvars -global
if obj.PreProcessingParameters.ROI_algorithm == 'PCA/ICA'
    ROI_tag = 1; % 1 for PCAICA, 2 for NMF
elseif obj.PreProcessingParameters.ROI_algorithm == 'NMF'
    ROI_tag = 2;
else
    error('Incorrect ROI algorithm')
end


% Now for PCAICA section
% first lets check if .icamasks file exists either due to prior clicking or
% due to running the CNN otherwise load .ica file
if ROI_tag == 1
    clearvars -global premasks
    poss_icamasks_name = strrep(obj.Dirs.runs{end}.sbx,'.sbx','.icamasks');
    if exist(poss_icamasks_name)
        load(poss_icamasks_name,'-mat')
        % will end up remaking the masks keep field so remove it
        premasks = rmfield(premasks,'masks_keep')
    else
        load(strrep(obj.Dirs.runs{end}.sbx,'.sbx','.ica'),'-mat');
    end
    % declare global variable and run GUI
    global premasks
    rr = ROI_selection_gui_PCAICA(premasks)
    uiwait(rr)
    % get output from gui and save
    global premasks
    if isfield(premasks,'masks_keep')
       save(poss_icamasks_name,'premasks','-v7.3')
    end
elseif ROI_tag == 2 % for CNMF
    clearvars -global premasks
    % load in the premasks variable
    % first check if already have a premasks variable clicked - via CNN for
    % example
    tmpR = dir([obj.Dirs.date_mouse '*.nmfmasks']);
    poss_premasks_filename = [obj.Dirs.date_mouse tmpR.name];
    if exist(poss_premasks_filename)
        load(poss_premasks_filename,'-mat')
        % will end up remaking the masks keep field so remove it
        premasks = rmfield(premasks,'masks_keep')
    else % haven't run with CNN or clicked before
        tmpR = dir([obj.Dirs.date_mouse '*.nmf']);
        poss_premasks_filename = [obj.Dirs.date_mouse tmpR.name];
        load(poss_premasks_filename,'-mat');
        min_ROI_size = 10;
        max_ROI_size = 500;
        premasks = PPPack.hf.throw_out_cells_too_small_too_large(premasks,min_ROI_size,max_ROI_size);
    end
    global premasks
    rr = ROI_selection_gui_NMF(premasks)
    uiwait(rr)
    % get output from gui and save
    global premasks
    if isfield(premasks,'masks_keep')
       save(poss_premasks_filename,'premasks','-v7.3')
    end
end

% update log file
obj.update_log_file('ROI_selection_GUI');
