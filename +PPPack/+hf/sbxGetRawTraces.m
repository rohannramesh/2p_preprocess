function [sig,signals] = sbxGetRawTraces(obj,dirs,premasks)

% for the size of the image and for appropriate information
z = PPPack.hf.sbxReadPMT(dirs.sbx,1,1);
sz = size(z);
% for the spatially downsampled size information
zz = PPPack.hf.sbxRemoveEdges(z);
siz_edge = size(zz);
zz = bin_mov_xyt(zz,2,1,0);
siz_ds = size(zz);
info = PPPack.hf.sbxInfo(dirs.sbx,1);

% make mask variable
signals = premasks.masks_keep;

% to protect against empty masks
for i = 1:length(signals)
    if isempty(signals(i).mask)
        signals(i).mask = zeros(size(signals(1).mask));
    end
end

% for the number of ROIs
nROIs = length(signals);


% bc clipped edges and spatially downsampled have to reverse this in the
% appropriate order
for curr_ROI = 1:nROIs
    curr_mask = signals(curr_ROI).mask;
    % first upsample by 2
    curr_mask = imresize(curr_mask,siz_edge,'nearest');
    % now add back edges
    curr_mask = PPPack.hf.sbxRemoveEdges([],curr_mask);
    signals(curr_ROI).mask = curr_mask;
end



if strcmp(obj.PreProcessingParameters.ROI_type, 'Axon') 
    signals = PPPack.hf.get_ROI_axons_masks_plus_neuropil(signals);
else
    signals = PPPack.hf.get_ROI_cellbody_masks_plus_neuropil(signals);
end
% update nROIs to include last ROI which is just neuropil
% rederive nROIs - and make variables that will replace as indexing
% approach
nROIs = length(signals);
mask = zeros(sz); % these are all binarized masks
mask_weights = zeros(sz); % these are all masks with weights
npil_masks = zeros(sz(1),sz(2),length(signals),'uint8'); % all neuropil masks
for curr_ROI = 1:nROIs
    mask = mask + curr_ROI.*(signals(curr_ROI).binmask);
    mask_weights = mask_weights + signals(curr_ROI).mask;
    npil_masks(:,:,curr_ROI) = uint8(signals(curr_ROI).neuropil);
end
% this variable will end up being the final variable - will save neuropil
% masks as a separate variable
signals_new.maskNumber = mask;
signals_new.maskWeights = mask_weights;


% pre-allocate for signal and neuropil
sig = zeros(info.max_idx+1, nROIs);
sig_npil = zeros(info.max_idx+1, nROIs);

use_registration_tag = 1;

%% Go multiple frames at a time with matlabpool and extract traces
% parallel processing
try a = gcp('nocreate');
catch err
    a = [];
end
if isempty(a)
    poolobj = parpool(10);
    pool_siz = poolobj.NumWorkers;;
else
    pool_siz = a.NumWorkers;
    poolobj = a;
end
nframes = info.max_idx+1;

chunk_size = 5000;
nchunks = ceil(nframes ./ chunk_size);
Frame_ind_start_points = [0:nchunks-1]*chunk_size;
nIter = ceil(nframes ./ (chunk_size*pool_siz));
RR_tmp_cell = zeros(nROIs,chunk_size,nchunks);
RR_tmp_npil = zeros(nROIs,chunk_size,nchunks);
if use_registration_tag == 1
    if strcmp(obj.PreProcessingParameters.registration,'whole pixel')
        ALIGNEDMAT = info.aligned.T;
    elseif strcmp(obj.PreProcessingParameters.registration,'subpixel')
        ALIGNEDMAT = info.aligned.StackReg;
    end
end
h = waitbar(0,sprintf('Straight pulling dem %d sigs',nIter));
fname = dirs.sbx;
for iBlock = 1:nIter
    %     waitbar(iBlock/(nIter),h);          % update waitbar...
    iFiles = intersect((iBlock-1)*pool_siz + [1:pool_siz], 1:nchunks)
    % if you have cleaned your data
    cleaned_data_path = [dirs.path dirs.date_mouse_run '_cleaned'];
    parfor curr_chunk_ind = iFiles
        if  ~exist(cleaned_data_path, 'dir')
            curr_chunk = Frame_ind_start_points(curr_chunk_ind);
            % if last chunk might not be 500 frame chunk
            if curr_chunk+chunk_size > nframes
                nframes_to_add_last_chunk = chunk_size - (curr_chunk+chunk_size-nframes);
                z = PPPack.hf.sbxReadPMT(fname,curr_chunk,nframes_to_add_last_chunk);
            else
                z = PPPack.hf.sbxReadPMT(fname,curr_chunk,chunk_size);
            end
            if use_registration_tag == 1
                if strcmp(obj.PreProcessingParameters.registration,'whole pixel')
                    for i = 1:size(z,3)
                        z(:,:,i) = circshift(z(:,:,i),ALIGNEDMAT(curr_chunk+i,:)); % align the image
                    end
                elseif strcmp(obj.PreProcessingParameters.registration,'subpixel')
                    if curr_chunk+chunk_size > nframes
                        [~,z] = PPPack.hf.sbxStackRegister(z,[],[],ALIGNEDMAT(curr_chunk+1:curr_chunk+nframes_to_add_last_chunk,:),0);
                    else
                        [~,z] = PPPack.hf.sbxStackRegister(z,[],[],ALIGNEDMAT(curr_chunk+1:curr_chunk+chunk_size,:),0);
                    end
                end
            end
            reshape_z = reshape(z,sz(1)*sz(2),size(z,3));
        else % if have denoised tiff data saved into cleaned folder
            mov = cleaned_data_path;
            list = dir(mov);             
            curr_tiff_name = list(curr_chunk_ind+2).name;
            [~, ~, ext] = fileparts(curr_tiff_name);
            if ~strcmp(ext, '.tif')
                continue
            end           
            z = PPPack.hf.readtiff(fullfile(mov, curr_tiff_name));
            reshape_z = reshape(z,size(z,1)*size(z,2),size(z,3));
        end
        for j=1:nROIs                          % for each cell and +1 for neuropil cell
            tmp_cell = find(signals_new.maskNumber == j); % this gives you the pixel values for that ROI
            tmp_cell_sig = mean(reshape_z(tmp_cell,:)); % this gets the mean for those pixels on that frame
            if size(tmp_cell_sig,2) ~= chunk_size % for last mov chunk that might not divide equally
                tmp_cell_sig(:,end+1:chunk_size) = 0; %put in +1
            end
            RR_tmp_cell(j,:,curr_chunk_ind) =  tmp_cell_sig;
            % neuropil
            tmp_npil = find(npil_masks(:,:,j));
            % if empty mask for neuropil - saying neuropil = 0 for the tc
            if isempty(tmp_npil)
                tmp_npil_sig = zeros(size(median(double(reshape_z(tmp_npil,:)))));
            else
                tmp_npil_sig = median(double(reshape_z(tmp_npil,:)));
            end
            if size(tmp_npil_sig,2) ~= chunk_size
                tmp_npil_sig(:,end+1:chunk_size) = 0; %put in +1
            end
            RR_tmp_npil(j,:,curr_chunk_ind) =  tmp_npil_sig;
        end
    end
end
delete(h)


if ~isempty(poolobj)
    delete(poolobj);
end

% now put back into easier format
for curr_ROI = 1:nROIs
    for curr_chunk_ind = 1:nchunks
        curr_chunk = Frame_ind_start_points(curr_chunk_ind);
        sig(curr_chunk+1:curr_chunk+chunk_size,curr_ROI) = RR_tmp_cell(curr_ROI,:,curr_chunk_ind);
        sig_npil(curr_chunk+1:curr_chunk+chunk_size,curr_ROI) = RR_tmp_npil(curr_ROI,:,curr_chunk_ind);
    end
end
% now make sure to clip off extra 0s at end (from matlabpool)
sig = sig(1:nframes,:);
sig_npil = sig_npil(1:nframes,:);
%% doing neuropil subtraction
sig_npil_subtracted = (sig - sig_npil);
median_sig = nanmedian(sig,1);
sig_npil_subtracted = bsxfun(@plus,sig_npil_subtracted,median_sig);

% Now put into signals format
for curr_ROI = 1:nROIs;
    signals_new.timecourse(curr_ROI).raw = sig(:,curr_ROI)';
    signals_new.timecourse(curr_ROI).neuropil = sig_npil(:,curr_ROI)';
    try if strcmp(obj.PreProcessingParameters.npil_correction,'Weighted')
            % get weight to scale npil to maximize skewness of subtracted trace
            subfun = @(x) -1 *skewness(signals_new.timecourse(curr_ROI).raw - ...
                (x * signals_new.timecourse(curr_ROI).neuropil));
            w = fminsearch(subfun,1);
            w(w<0) = 0; % can't be less than 0 or greater than 2
            w(w>1.5) = 1;
%             w(w>2) = 2; % used to be 2 and used to just add back median
%             of raw - but now adding back weight * npil
            signals_new.timecourse(curr_ROI).subtracted = (signals_new.timecourse(curr_ROI).raw  - ...
                (w.*signals_new.timecourse(curr_ROI).neuropil))+w*nanmedian(signals_new.timecourse(curr_ROI).neuropil);
            signals_new.npil_weight(curr_ROI) = w;
        else
            signals_new.timecourse(curr_ROI).subtracted = sig_npil_subtracted(:,curr_ROI)';
            signals_new.npil_weight(curr_ROI) = 1;
        end
    catch err
        signals_new.timecourse(curr_ROI).subtracted = sig_npil_subtracted(:,curr_ROI)';
        signals_new.npil_weight(curr_ROI) = 1;
    end
end

% now lets redefine the signals variable
clear signals
signals = signals_new;

% save the npil masks
save(strrep(dirs.sbx,'.sbx','.npilmasks'),'npil_masks','-v7.3')

% save([dirs.sbx '.signals'],'sig','signals');     % append the motion estimate data...

