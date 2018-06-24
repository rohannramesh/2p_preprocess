function obj = sbxRegister(obj)
% register all runs from a given directory using efficient subpixel
% registration (Bonin et al. (2011)). The output of this file will be 
% a .align file saved into each folder. This will contain the
% subpixel registration output called StackReg. 
% Columns for this variable:
%      OUTS(:,1) is correlation coefficients
%      OUTS(:,2) is global phase difference between images (target and
%                image to be registered
%                (should be zero if images real and non-negative).
%      OUTS(:,3) is net row shift
%      OUTS(:,4) is net column shift
% If you need to apply the registration use the function PPPack.hf.sbxStackRegister
% However, to load aligned data, I recommend you use sbxLoadReg or
% sbxLoadRegRun within PPPack.hf

% iterate through each run
nRuns = length(obj.Dirs.runs);
for curr_run = 1:nRuns
    % this is the directory information for a specific run
    dirs = obj.Dirs.runs{curr_run};

    % quickly load a single frame for size and sbxInfo purposes
    z = PPPack.hf.sbxReadPMT(dirs.sbx, 1, 1);

    % if already aligned this run then don't continue registering the data
    if exist(strrep(dirs.sbx,'.sbx','.align'))
        continue
    end

    % This is the info that is stored in the default sbx file that tells
    % you the number of frames, the framerate, the size of the image etc
    info = PPPack.hf.sbxInfo(dirs.sbx,1);

    % define number of frames
    nframes = info.max_idx;

    % settings to use for subpixel registration - DON'T TOUCH
    defaultopts = {'Overwrite',true,'Oversampling',4,'DoRecurse',true,...
                    'Engine','subpixel','TargetFileNs',1,'PoolSize',1,'recgreenpn_ref',[]};
    options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);
    options
    os = options.Oversampling;

    % target is 500 frames of activity from 200 - 700 from the target run
    % register this once internally and then take the mean as the target
    % for all future runs
    first_frame = 200;
    if isempty(dirs.target) 
        display(['Registering to ' dirs.sbx])
        target_stack = PPPack.hf.sbxReadPMT(dirs.sbx, first_frame-1, 500);
    else
        display(['Registering to ' dirs.target])
        target_stack = PPPack.hf.sbxReadPMT(dirs.target,first_frame-1,500); 
    end
    target_backup2 = target_stack;
    % register only internal part of image bc of blacked out edges from the
    % imaging rig itself
    frac1 = .15; % Fraction of image to keep - used to be 0.25 - only the internal portion of the image
    sz0 = size(target_backup2);
    target_stack_small = target_backup2(floor(sz0(1)*frac1):ceil(sz0(1)*(1-frac1)), ...
        floor(sz0(2)*frac1):ceil(sz0(2)*(1-frac1)),:);
    % this is the function that does the registration
    [~,target_stack_reg] = PPPack.hf.sbxStackRegister(target_stack_small,mean(target_stack_small,3),os);
    target = mean(target_stack_reg,3); % this is the single mean image all runs will be registered to
    target_backup = target;



    %% for loop with matlabpool to parallelize registration
    a = gcp('nocreate');
    if isempty(a)
        poolobj = parpool(10);
        pool_siz = poolobj.NumWorkers;
    else
        pool_siz = a.NumWorkers;
        poolobj = a;
    end
    
    % We will divide our run into 500 frame chunks so that we can
    % parallelize processing using parfor loops
    chunk_size = 500;
    nchunks = ceil(nframes ./ chunk_size);
    Frame_ind_start_points = [0:nchunks-1]*chunk_size; % these are the frames to start each chunk
    nIter = ceil(nframes ./ (chunk_size*pool_siz));
    ChunkVar = zeros(500,4,nchunks); % this will be where we save the shifts to use for registration
    T = zeros(nframes,2);
    for iBlock = 1:nIter
        iFiles = intersect((iBlock-1)*pool_siz + [1:pool_siz], 1:nchunks)
        % parfor loop to load in each chunk of data register it and then
        % save the shifts necessary for saving
        parfor curr_chunk_ind = iFiles
            curr_chunk = Frame_ind_start_points(curr_chunk_ind);
            % if last chunk might not be 500 frame chunk
            % this is the section actually loading in the movie in the
            % appropriate chunks
            if curr_chunk+chunk_size > nframes
                nframes_to_add_last_chunk = chunk_size - (curr_chunk+chunk_size-nframes)+1;
                z = PPPack.hf.sbxReadPMT(dirs.sbx,curr_chunk,nframes_to_add_last_chunk);
            else
                z = PPPack.hf.sbxReadPMT(dirs.sbx,curr_chunk,chunk_size);
            end
            % we crop the chunk of movie to match the cropped target
            z = z(floor(sz0(1)*frac1):ceil(sz0(1)*(1-frac1)), ...
                floor(sz0(2)*frac1):ceil(sz0(2)*(1-frac1)),:);
            try [shifts] = PPPack.hf.sbxStackRegister(z,target,os);
                if size(shifts,1) ~= size(ChunkVar(:,:,curr_chunk_ind),1)
                    shifts(end+1:size(ChunkVar(:,:,curr_chunk_ind),1),:) = 0; 
                end
                ChunkVar(:,:,curr_chunk_ind) = shifts;
            catch err
            end
        end
    end
    T = [];
    for i = 1:nchunks
        T = [T;ChunkVar(:,:,i)];
    end
    % remove any of padded frames at end
    ind_toss = find(T(:,1) == 0 & T(:,2) == 0);
    T(ind_toss,:) = [];
    StackReg = T;
    T = round(T(:,3:4)); % needs to be an integer number for their reg using circshift


    if ~isempty(poolobj)
        delete(poolobj);
    end
    %% Now save data
    m = uint16(target_backup); % mean image of target
    if isfield(dirs, 'target')
        registered_to = dirs.target;
    else
        registered_to = [];
    end
    save([dirs.path '\' dirs.sbx_name '.align'],'m','T','StackReg','registered_to');
end

obj.update_log_file('registration')
