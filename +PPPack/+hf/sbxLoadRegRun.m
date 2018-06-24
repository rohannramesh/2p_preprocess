function mov = sbxLoadRegRun(dirs,reg_type)
% Building in the -1 bc of 0 indexing - assuming the frame you give this
% function is not zero indexed - so will load in the frame 1 back in sbx
% language


% set curr dir
cd(dirs.path)

% get filename
fname = dirs.sbx(1:end-4);

% load .align file
RegMat = load([fname '.align'],'-mat');

info = PPPack.hf.sbxInfo(dirs.sbx);

% if have a reg sbx file already written
if ~isempty(dirs.sbxreg)
    display('Found .sbxreg file and using that')
    mov = PPPack.hf.sbxReadPMT(dirs.sbxreg,0,info.max_idx+1);
else

    % define chunking parameters
    chunk_size = 5000;
    chunks = 1:chunk_size:info.max_idx+1; corenum = 10;
    nchunks = numel(chunks); % number of indexes
    chunk_idx = arrayfun(@(i) chunks(i:min((i+corenum-1), nchunks)), ...
        1:corenum:nchunks, 'UniformOutput', false); % indices of the patches in each batch
    % pre allocate
    mov_chunks = cell(1,nchunks);
    
    parfor ii = 1:nchunks
        idx = chunks(ii);

        % load frames
        nframes = size(RegMat.T,1);
        idx:min(idx+chunk_size-1,nframes);
        if idx+chunk_size > nframes
            mov_chunks{ii} = PPPack.hf.sbxReadPMT(dirs.sbx,idx-1,nframes - idx + 1);
        else
            mov_chunks{ii} = PPPack.hf.sbxReadPMT(dirs.sbx,idx-1,chunk_size);
        end
        % squeeze to get rid of other PMTs
        mov_chunks{ii} = squeeze(mov_chunks{ii});
        % get info variabl

        nframes = info.max_idx+1;
        % get correct registration technique
        if strcmp(reg_type,'whole pixel')
            ALIGNEDMAT = info.aligned.T;
        elseif strcmp(reg_type,'subpixel')
            ALIGNEDMAT = info.aligned.StackReg;
        end
        if strcmp(reg_type,'whole pixel')
            for i = 1:size(mov_chunks{ii},3)
                mov_chunks{ii}(:,:,i) = circshift(mov_chunks{ii}(:,:,i),ALIGNEDMAT(idx+i-1,:)); % align the image
            end
        elseif strcmp(reg_type,'subpixel')
            if idx+chunk_size > nframes
                [~,mov_chunks{ii}] = PPPack.hf.sbxStackRegister(mov_chunks{ii},[],[],ALIGNEDMAT(idx:end,:),0);
            else
                [~,mov_chunks{ii}] = PPPack.hf.sbxStackRegister(mov_chunks{ii},[],[],ALIGNEDMAT(idx:idx+chunk_size-1,:),0);
            end
        end
    end

    % unshuffle
    mov = zeros(info.sz(1),info.sz(2),info.max_idx+1,'uint16');
    ind_full = 1;
    for i = 1:length(mov_chunks)
        mov(:,:,ind_full:min(ind_full+chunk_size-1,info.max_idx+1)) = mov_chunks{i};
        mov_chunks{i} = [];
        ind_full = ind_full + chunk_size;
    end
    
end


