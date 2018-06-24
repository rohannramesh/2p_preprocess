function obj = ROI_selection_CNN_NMF(obj)
% Select which of the putative ROIs identified by the NMF algorithm
% to keep for ROI trace extraction via use of a pre-trained Convolutional
% Neural Network. This neural network was pre-trained by Rohan Ramesh from
% over 50 imaging sessions of manually hand-selected ROIs.
% this function will create a folder with a images in them for reading into
% the CNN. The output of the CNN will be a set of probabilities
% corresponding to whether the putative ROI should be considered a neuron
% [1 0] or not [0 1]
% In order to run must have installed Python 3.5, Numpy, Keras, and
% tensorflow and loaded into matlab via the command pyversion
% The output will be a Predictions.npy variable and builds a .nmfmasks
% variable saved in the day directory (e.g. 170524_RR) that contains those
% masks to keep according to the CNN

if isempty(pyversion)
    error('Need to register python with matlab using pyversion')
end

% threshold for calling a neuron a neuron and minimum/maximimum size for an ROI 
CNN_threshold = 0.5;
min_ROI_size = 10;
max_ROI_size = 500;

% Necessary paths to know where the pretrained CNN and the run python
% script is
p = mfilename('fullpath');
ind_slash = find(p == '\');
basepath_to_use = p(1:ind_slash(end-1));
% CNN_path = 'D:\Analysis_scripts\Dropbox\AndermannLab\users\ramesh\+PPPack\+PCAICA\PCAICA_CNN.h5';
CNN_path = [basepath_to_use '+NMF\NMF_CNN.h5'];
% python script path
% python_script_CNN = 'D:\Analysis_scripts\Dropbox\AndermannLab\users\ramesh\+PPPack\+PCAICA\run_CNN_python.py';
python_script_CNN = [basepath_to_use '+NMF\run_CNN_NMF_python.py'];

% first need to load in premasks variable
tmpR = dir([obj.Dirs.date_mouse '*.nmf']);
premasks_filename = [obj.Dirs.date_mouse tmpR.name];
load(premasks_filename,'-mat');
% % resave with legacy version
% save(strrep(premasks_filename,'.nmf','.legacynmf'),'premasks',...
%     'A','C','keep','b','f','S','P','RESULTS','YrA','ROIvars','-v7.3')
% lets rederive premasks but with limits in terms of number of pixel size
premasks = PPPack.hf.throw_out_cells_too_small_too_large(premasks,min_ROI_size,max_ROI_size);

% make tiff folder with all tiffs for CNN testing
savefolder_tiffs = [obj.Dirs.date_mouse,'ROI_images'];
PPPack.hf.make_tiff_folder_ROIs_CNN_multiple_erosions(premasks,savefolder_tiffs);
% this is the path to the .mat file with the histograms
tc_hist_path = [savefolder_tiffs '\tc_hist.mat'];

% running
% had to run using system and powershell to deal with python crashes bc of
% low-level administrative write problems - built for windows server system
system(['powershell -Command Start-Process python.exe '...
    '-ArgumentList "' python_script_CNN '","' CNN_path '","' savefolder_tiffs '","' tc_hist_path '" -Verb RunAs'])
pause(40) % bc takes some time to run

% lets read the .npy file and choose which ROIs to keep
CNN_output = PPPack.hf.readNPY([savefolder_tiffs '\Predictions.npy']);
% get all of the ROI names and erosion levels from that folder 
tif_names = dir([savefolder_tiffs '\All\*.tif']);
if length(tif_names) ~= size(CNN_output,1)
    error('Different number of tiffs from CNN output')
end

% build variable that includes ROI number the erosion level and the CNN
% probability bc this will allow us to resort the ROIs based on the
% threshold for 
ROIs_erosion_CNN = [];
for i = 1:size(CNN_output,1)
    ind_und = find(tif_names(i).name == '_');
    ind_dot = find(tif_names(i).name == '.');
    ROI_number = str2num(tif_names(i).name(ind_und(1)+1:ind_und(2)-1));
    erosion_level = str2num(tif_names(i).name(ind_und(end)+1:ind_dot-1))./10; % divide by 10 bc scaled up in tiff making
    ROIs_erosion_CNN = [ROIs_erosion_CNN; ROI_number erosion_level CNN_output(i,1)];    
end
% lets take the ROI and the erosion level that gave the best CNN result
unique_ROIs_to_keep = unique(ROIs_erosion_CNN(:,1));
for curr_ROI = unique_ROIs_to_keep'
    idx = find(ROIs_erosion_CNN(:,1) == curr_ROI);
    best_one = find(ROIs_erosion_CNN(idx,3) == max(ROIs_erosion_CNN(idx,3)));
    if length(best_one) > 1
        best_one = best_one(end); % if same CNN probability take biggest mask
    end
    ROIs_erosion_CNN(setdiff(idx,idx(best_one)),:) = [];
end
% resort by the CNN threshold
[~,CNN_order] = sort(ROIs_erosion_CNN(:,3),'descend');
ROIs_erosion_CNN = ROIs_erosion_CNN(CNN_order,:);
% now resave the ica file appropriately sorted
premasks.ica_backup = premasks.ica;
premasks = rmfield(premasks,'ica'); % remove before overwriting
for curr_ROI_idx = 1:size(ROIs_erosion_CNN)
    sorted_ROI = ROIs_erosion_CNN(curr_ROI_idx,1);    
    premasks.ica(curr_ROI_idx) = premasks.ica_backup(sorted_ROI);
    % update so using the eroded mask but only if supposed to keep it
    if ROIs_erosion_CNN(curr_ROI_idx,3) > CNN_threshold
        premasks.ica(curr_ROI_idx).filter = PPPack.hf.get_erosion_mask(...
            premasks.ica(curr_ROI_idx).filter,ROIs_erosion_CNN(curr_ROI_idx,2));
    else
         premasks.ica(curr_ROI_idx).filter = premasks.ica(curr_ROI_idx).filter;
    end
    % update the idx variable
    premasks.ica(curr_ROI_idx).idx = [];
    [premasks.ica(curr_ROI_idx).idx(:,1) premasks.ica(curr_ROI_idx).idx(:,2)] = ...
        find(single(premasks.ica(curr_ROI_idx).filter>0));
end
new_ROI_idx_to_keep = find(ROIs_erosion_CNN(:,3) > CNN_threshold);
% reset ROI list bc using CNN
premasks.ROI_list = zeros(size(premasks.ROI_list));
premasks.ROI_list(new_ROI_idx_to_keep) = 1;
premasks.erosions_used_CNN = ROIs_erosion_CNN(:,2);
% remove backup field
premasks = rmfield(premasks,'ica_backup');
% just a heads up that did CNN
premasks.did_CNN = 1;
% bc did CNN and NMF already populates these fields have to reset mask image fields
premasks.masks = zeros(size(premasks.masks));
premasks.masks_with_group = zeros(size(premasks.masks));




% lets edit the premasks variable with the suggestions from the CNN and
% create an output for either straight extraction of traces or of viewing
% within the GUI
premasks = PPPack.hf.build_masks_for_extraction(premasks);

% remove doubles - ROIs that would be double counted
premasks = PPPack.hf.remove_double_masks(premasks);


% now add the masks to appropriate variables
ROI_idx_from_total_ica = find(premasks.ROI_list);
for i = 1:length(premasks.masks_keep)
    curr_ROI = i;
    premasks.masks = premasks.masks + ...
        ROI_idx_from_total_ica(curr_ROI)*single(premasks.masks_keep(curr_ROI).mask>0);
    % masks with group number
    premasks.masks_with_group = premasks.masks_with_group + ...
        1*single(premasks.masks_keep(curr_ROI).mask>0);
    % can't have over 5 groups
    premasks.masks_with_group(premasks.masks_with_group > 5) = 5;
end

% remove the directory containing the tiffs but first move the predictions
% file
movefile([savefolder_tiffs '\Predictions.npy'],[obj.Dirs.date_mouse 'Predictions.npy']);
rmdir(savefolder_tiffs,'s')


% save this variable
save(strrep(premasks_filename,'.nmf','.nmfmasks'),'premasks','-v7.3')

% update log file
obj.update_log_file('CNN_NMF');