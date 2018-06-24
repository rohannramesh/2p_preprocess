function obj = CNMF_ROI_identification(obj)
% Run Constrained Non-Negative Matrix Factorization ROI
% identification algorithm to find putative cell masks
% This algorithm will divide up the movie into spatial and temporal
% matrices (A and C respectively that correspond to ROIs and
% deconvolved timecourses (implemented by Rohan Ramesh, see
% Pnevmatikakis et al., 2016) 
% Output of running this will save a *.nmf file in the directory of the
% day folder. This file will be formatted to either be run through
% the Matlab GUI, preprocessed for the Web Clicking GUI, or run through a
% pre-trained Convolutional Neural Network for ROI selection


nRuns = length(obj.Dirs.runs);
% create name for matfile saving purposes - will save this in the folder
% outside the last run to prevent confusion with PCAICA analysis
data_matfile_name = ['data_runs'];
runs_appendage = ['runs'];
for curr_run = 1:nRuns
    tmpDir = obj.Dirs.runs{curr_run};
    run_number_idx = strfind(tmpDir.date_mouse_run,'run');
    run_number = tmpDir.date_mouse_run(run_number_idx+3:end);
    data_matfile_name = [data_matfile_name num2str(run_number)];
    runs_appendage = [runs_appendage num2str(run_number)];
end

% first check to see if already built .matfile and if so don't build again
data_matfile_path = sprintf('%s\%s.mat',obj.Dirs.date_mouse,data_matfile_name);
if ~exist(data_matfile_path)
    mov_ds = [];
    for curr_run = 1:nRuns
        % this is the directory information for a specific run
        curr_dir = obj.Dirs.runs{curr_run};

        % for now also cd into the proper directory
        cd(curr_dir.path)
        % load in ds mov
        mov =  PPPack.hf.sbxLoadDSMov(obj,curr_dir);

        sz_mov = size(mov);
        % concatenate
        mov_ds = cat(3,mov_ds,mov);
    end

    % now lets set up matfile
    data = matfile(data_matfile_path,'Writable',true);
    % lets put the variables into it
    data.Y = mov_ds; % the entire movie
    sizY = [size(mov_ds)]; % size of the movie
    data.Yr = reshape(mov_ds,prod(sizY(1:end-1)),[]); % each column is single frame linearized
    data.nY = min(mov_ds(:));
    data.sizY = sizY;
    data.movm = mean(mov_ds,3);
    mov_tds = PPPack.hf.group_z_project(mov_ds,100);
    data.CrossCorrMeanImg = PPPack.hf.CrossCorrImage(mov_tds);
    clear mov_ds mov_tds
else
    data = matfile(data_matfile_path,'Writable',true);
end

%% Run CNMF algorithm - should be on your path already if you initialized your package correctly and then added it to path

% Set parameters
sizY = size(data,'Y');                                  % size of data matrix
patch_size = obj.PreProcessingParameters.Patch_Size;    % size of each patch along each dimension (optional, default: [32,32])
% patch_size = [60,60];                                 % size of each patch along each dimension (optional, default: [32,32])
overlap = [8,8];                                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 40;                                           % number of components to be found
tau = obj.PreProcessingParameters.Cell_Width;     % std of gaussian kernel (size of neuron) - used to be 8 - 5
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold
sizY = data.sizY;

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'min_corr',0.1,...
    'search_method','ellipse','dist',2,...      % search locations when updating spatial components - used to be 3 RR
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'cluster_pixels',false,...
    'ssub',1,...
    'tsub',1,...
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',tau,... 
    'spatial_method','constrained',...
    'min_size',1,...
    'gSig',1,...
    'min_pixel',3,...
    'min_size_thr',3,...
    'max_timesteps',sizY(3));

%% Run on patches

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

%% classify components
[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(data,A,C,b,f,YrA,options);

%% add tc to ROIvars

ROIvars = PPPack.hf.add_tc_to_ROIvars(data,A,b,ROIvars);

%%
global premasks
[premasks] = PPPack.hf.convert_nmf_to_premasks(obj.Dirs,data,A,C,keep,ROIvars)

%%
[curr_mouse,curr_date] = PPPack.hf.get_mouse_day_run_info_from_dirs(curr_dir);
save([obj.Dirs.date_mouse curr_mouse '_' curr_date '_' runs_appendage '.nmf'],'premasks',...
    'A','C','keep','b','f','S','P','RESULTS','YrA','ROIvars','-v7.3')

% update log file
obj.update_log_file('NMF_ROI_identification');
