function obj = sbxTestRegistrationImageStack(obj)
% Generates tif file from sbx files
% registers images if .align file has been created
% saves tif in 500 frame chunks in new folder if obj.PreProcessingParameters.write_first_last_500 = 0
% if obj.PreProcessingParameters.write_first_last_500 = 1 then only writes first and last 500 frames

% iterate through each run
nRuns = length(obj.Dirs.runs);
for curr_run = 1:nRuns
    % this is the directory information for a specific run
    dirs = obj.Dirs.runs{curr_run};

    % this is to get the size of the image
    fname = dirs.sbx;
    z = PPPack.hf.sbxReadPMT(fname, 1, 1);
    sz = size(z(1,:,:));
    % to get the info from a given run - IF already aligned this will have
    % the registration shifts necessary for registration of movie
    info = PPPack.hf.sbxInfo(dirs.sbx,1);

    % this tag makes sure you only write the first and last X number of
    % frames to test registration as opposed to writing the entire tiff
    % stack
    write_only_first_last_tag = obj.PreProcessingParameters.write_first_last_500;


    N = info.max_idx;




    %% write first and last mov frame
    chunk_size = 500;
    nchunks_need_to_write = ceil(N./chunk_size);
    if write_only_first_last_tag == 1
        write_first_last_chunk(fname,chunk_size,info,N,obj);
        return
    end

    % make new directory to save tif
    if ~isempty(info.aligned)
        tif_dir = [dirs.sbx(1:end-4) '_tiff_stack_reg'];
    else
        tif_dir = [dirs.sbx(1:end-4) '_tiff_stack'];
    end
    if ~exist(tif_dir)
        mkdir(tif_dir);
    end

    %% write whole movie
    k = 0;
    done = 0;
    chunk_ind = 1;
    new_tif = 0;
    curr_stack_ind = 1;

    % write first and last 500
    if write_only_first_last_tag == 1
        write_first_last_chunk(fname,chunk_size,info,N,obj);
    end
    % now write whole tiff stack
    while curr_stack_ind <= chunk_size && k <=N
        if new_tif == 1
            curr_stack = nan(size(z,1),sz(2),chunk_size); 
        end    
        q = PPPack.hf.sbxReadPMT(fname, k, 1);
        % just write first channel
        % shift by x and y either using subpixel or whole pixel shifts
        if ~isempty(info.aligned)
            if strcmp(obj.PreProcessingParameters.registration,'whole pixel')
                q = circshift(q,info.aligned.T(k+1,:)); % align the image
            elseif strcmp(obj.PreProcessingParameters.registration,'subpixel')
                [~,q] = PPPack.hf.sbxStackRegister(q,[],[],info.aligned.StackReg(k+1,:),0);
            else 
                error('Unidentified registration type')
            end
        end
        curr_stack(:,:,curr_stack_ind) = q;
        curr_stack_ind = curr_stack_ind + 1;
        new_tif = 0;
        k = k + 1;
        if curr_stack_ind > chunk_size
            new_tif = 1;
            curr_stack_ind = 1;
            % now write tif
            PPPack.hf.writetiff(curr_stack,[tif_dir '\' fname '_' num2str(chunk_ind,'%04.0f') '.tif'],'uint16');
            chunk_ind = chunk_ind + 1;
        end
    end
end
end

function write_first_last_chunk(fname,chunk_size,info,N,obj)
    % load in first chunk of frames
    ff = PPPack.hf.sbxReadPMT(fname, 0, chunk_size);
    if ~isempty(info.aligned)
        if strcmp(obj.PreProcessingParameters.registration,'whole pixel')
            for i = 1:size(ff,3)
                ff(:,:,i) = circshift(ff(:,:,i),info.aligned.T(1+i-1,:)); % align the image
            end
        elseif strcmp(obj.PreProcessingParameters.registration,'subpixel')   
            [~,ff] = PPPack.hf.sbxStackRegister(ff,[],[],info.aligned.StackReg(1:chunk_size,:),0);
        end
    end
    % load in last chunk of frames
    lf = PPPack.hf.sbxReadPMT(fname, N-chunk_size+1, chunk_size);
    if ~isempty(info.aligned)
        if strcmp(obj.PreProcessingParameters.registration,'whole pixel')
            for i = 1:size(lf,3)
                lf(:,:,i) = circshift(lf(:,:,i),info.aligned.T(N-chunk_size+1+i,:)); % align the image
            end
        elseif strcmp(obj.PreProcessingParameters.registration,'subpixel')   
            [~,lf] = PPPack.hf.sbxStackRegister(lf,[],[],info.aligned.StackReg(N-chunk_size+2:N+1,:),0);
        end
    end
    total = cat(3,ff,lf);
    PPPack.hf.writetiff(total,[fname '_first_last_' num2str(chunk_size) '_frames'],'uint16')
end

