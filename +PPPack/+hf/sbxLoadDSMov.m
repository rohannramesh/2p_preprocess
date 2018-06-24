function mov = sbxLoadDSMov(obj,dirs)
% this will be the value we downsample by
ds_val = 5; % in time
spatial_ds = 2;
% we are going to spatially downsample by a factor of 2 to make things
% faster - but first need to clip the edges of the movie
c = PPPack.hf.sbxReadPMT(dirs.sbx, 1, 1);
c = PPPack.hf.sbxRemoveEdges(c);
c = PPPack.hf.bin_mov_xyt(c,spatial_ds,1,0);
sz = size(c);
global info;

% number of frames
N = info.max_idx+1;


if ~isempty(info.aligned)
    display(['Giving registered movie and using ' obj.PreProcessingParameters.registration ' registration'])
    use_registration_tag = 1;
else
    use_registration_tag = 0;
end

%% Matlabpool version

%% for loop with matlabpool
a = gcp('nocreate');
if isempty(a)
    poolobj = parpool(10);
    pool_siz = poolobj.NumWorkers;
else
    pool_siz = a.NumWorkers;
    poolobj = a;
end

% decide chunk size for loading in and registering + crop, spatially
% downsample, and then temporally downsample
chunk_size = 5000;
eff_chunk_size = chunk_size./ds_val;
% pre allocation
mod_val = ceil(N/chunk_size)*chunk_size;
mov = zeros(sz(1),sz(2),mod_val./ds_val);
nframes = N;
nchunks = ceil(nframes ./ chunk_size);
Frame_ind_start_points = [0:nchunks-1]*chunk_size;
nIter = ceil(nframes ./ (chunk_size*pool_siz));
RR_tmp_cell = zeros(sz(1),sz(2),chunk_size./ds_val,nchunks);
T = zeros(nframes,2);
if use_registration_tag == 1
    if strcmp(obj.PreProcessingParameters.registration,'whole pixel')
        ALIGNEDMAT = info.aligned.T;
    elseif strcmp(obj.PreProcessingParameters.registration,'subpixel')
        ALIGNEDMAT = info.aligned.StackReg;
    end
end
% iterate through load movie in, crop, spatially downample, and then
% temporally downsample
for iBlock = 1:nIter
%     waitbar(iBlock/(nIter),h);          % update waitbar...
    iFiles = intersect((iBlock-1)*pool_siz + [1:pool_siz], 1:nchunks)
    parfor curr_chunk_ind = iFiles
        curr_chunk = Frame_ind_start_points(curr_chunk_ind);
        % if last chunk might not be 5000 frame chunk
        if curr_chunk+chunk_size > nframes
            nframes_to_add_last_chunk = chunk_size - (curr_chunk+chunk_size-nframes);
            z = PPPack.hf.sbxReadPMT(dirs.sbx,curr_chunk,nframes_to_add_last_chunk);
        else
            z = PPPack.hf.sbxReadPMT(dirs.sbx,curr_chunk,chunk_size);
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
        % first lets crop
        z = PPPack.hf.sbxRemoveEdges(z);        
        % spatially downsample by factor of 2
        z = PPPack.hf.bin_mov_xyt(z,spatial_ds,1,0);
        if size(z,3) ~= chunk_size
            z(:,:,end+1:chunk_size) = 0; 
        end
        % downsample temporally
        z = PPPack.hf.downsampleWithAvg(z,ds_val);
        RR_tmp_cell(:,:,:,curr_chunk_ind) = z;
    end
end

% now put into standard mov format
ind = 1;
for curr_chunk_ind = 1:nchunks
    mov(:,:,ind:ind+eff_chunk_size-1) = (RR_tmp_cell(:,:,:,curr_chunk_ind));
    ind = ind+eff_chunk_size;
end
mov = mov(:,:,1:(N./ds_val));


if ~isempty(poolobj)
    delete(poolobj);
end




