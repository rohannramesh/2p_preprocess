% demo script for splitting the field of view in patches and processing in parallel
% through memory mapping. See also run_pipeline.m for the complete  
% pre-processing pipeline of large datasets

clear;
%% load file

path_to_package = '../ca_source_extraction';   % path to the folder that contains the package
addpath(genpath(path_to_package));
             
filename = '\\megatron\Analysis_scripts\F1392_slice2_teststack_rawF.tif';      % path to stack tiff file
%foldername = '';    % path to folder that contains a sequence of tiff files

% if exist([filename(1:end-4),'v2','.mat'],'file')
% %     data = matfile([filename(1:end-3),'mat'],'Writable',true);
%     data = matfile([filename(1:end-4),'v2','.mat'],'Writable',true);
% else
%     sframe=1;						% user input: first frame to read (optional, default 1)
%     num2read=[];					% user input: how many frames to read   (optional, default until the end)
% %     chunksize=5000;                 % user input: read and map input in chunks (optional, default read all at once)
%     chunksize=[];                 % user input: read and map input in chunks (optional, default read all at once)
%     data = memmap_file(filename,sframe,num2read,chunksize);
% %     data = memmap_file_sbx(filename,sframe,num2read,chunksize);
%     %data = memmap_file_sequence(foldername);
% end
    Y = load_tiff(filename);
    sizY = size(Y);
    Yr = reshape(Y,prod(sizY(1:end-1)),[]);
    nY = min(Yr(:));
    %Yr = Yr - nY;
    %save([filename(1:end-3),'mat'],'Yr','Y','nY','sizY','-v7.3');
%     savefast([filename(1:end-3),'mat'],'Yr','Y','nY','sizY');
    data = matfile([filename(1:end-3),'mat'],'Writable',true);
%     data = matfile([filename(1:end-3),'mat'],'Writable',true);
    data.Y = Y;
    data.Yr = Yr;
    data.nY = nY;
    data.sizY = [sizY(1:end-1),size(Y,3)-sframe+1];    
    data.movm = mean(Y,3);
    mov_tds = group_z_project(Y,100);
    data.CrossCorrMeanImg = CrossCorrImage(mov_tds);
    clear mov_ds mov_tds
%% Set parameters
%% Set parameters
sizY = size(data,'Y');                    % size of data matrix
% patch_size = params.Patch_Size;           % size of each patch along each dimension (optional, default: [32,32])
patch_size = [30,30];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [8,8];                          % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 40;                                            % number of components to be found
tau = 3;                            % std of gaussian kernel (size of neuron) - used to be 8 - 5
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold
sizY = data.sizY;

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'min_corr',0.1,...
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
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

ROIvars = add_tc_to_ROIvars(data,A,b,ROIvars);
%%
% change for save location
MultiRunDirs.date_mouse = '\\megatron\Analysis_scripts'
curr_mouse = 'TEST'
curr_date = 'TEST2'

global icaguidata
[icaguidata] = convert_nmf_to_icaguidata(MultiRunDirs,data,A,C,keep,ROIvars)

%%
% [curr_mouse,curr_date] = get_mouse_date_runs_from_path(MultiRunDirs.runs{1}.path);
save([MultiRunDirs.date_mouse '\' curr_mouse '_' curr_date '.icanmf'],'icaguidata',...
    'A','C','keep','b','f','S','P','RESULTS','YrA','ROIvars','-v7.3')

%% cell-clicking

rr = cellSortChooseRoiGUI_nmf(icaguidata)

    
%%
% now re-extract the icaguidata file and save it
global icaguidata
% save wherever you want NOW
