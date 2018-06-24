function make_tiff_folder_ROIs_CNN_multiple_erosions(premasks,directory_save)
directory_save_root = directory_save;
directory_save = [directory_save '\All'];

% first lets make the directory if doesnt exist
if ~exist(directory_save)
    mkdir(directory_save)
end
nROIs = length(premasks.ica);
minim = 20;
if isfield(premasks,'ROIs_selected_from_orig') % this means from NMF
    all_erosion_levels = [1 0.9 0.8 0.7];
else
    all_erosion_levels = [1 0.7 0.5 0.3 0.1];
end
tc_bins_use = 0:0.01:1;
all_hists = [];
for curr_ROI = 1:nROIs
    curr_mask = premasks.ica(curr_ROI).filter;
    % make images
    for erosion_level = all_erosion_levels
        if erosion_level ~= 1
            new_mask = PPPack.hf.get_erosion_mask(curr_mask, erosion_level);
        else
            new_mask = curr_mask;
%             display('KELLY IS WRONG')
        end
        if sum(sum(new_mask)) == 0
            continue
        end
        impath = [directory_save '\filter_' num2str(curr_ROI, '%05d') '_er_' num2str(erosion_level*10,'%03d') '.tif'];
        PPPack.hf.saveFilterImage_forCNN(impath, new_mask, minim); 
        % make timecourse histogram
        curr_hist = PPPack.hf.get_hist_of_tc(premasks,curr_ROI,tc_bins_use);
        all_hists = [all_hists; curr_hist];
    end

end
% now save all hists
save([directory_save_root '\tc_hist'],'all_hists')
end


