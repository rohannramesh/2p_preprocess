function make_tiff_folder_ROIs_CNN(premasks,directory_save)

directory_save = [directory_save '\All'];

% first lets make the directory if doesnt exist
if ~exist(directory_save)
    mkdir(directory_save)
end
nROIs = length(premasks.ica);
minim = 20;
for curr_ROI = 1:nROIs
    impath = [directory_save '\filter_' num2str(curr_ROI, '%05d') '.tif'];
    PPPack.hf.saveFilterImage_forCNN(impath, premasks.ica(curr_ROI).filter, minim); 
end
