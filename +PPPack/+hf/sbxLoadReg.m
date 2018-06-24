function mov = sbxLoadReg(dirs,reg_type,start_frame,chunk_size)
% Building in the -1 bc of 0 indexing - assuming the frame you give this
% function is not zero indexed - so will load in the frame 1 back in sbx
% language


% % set curr dir
% cd(dirs.path)

% load .align file
RegMat = load(strrep(dirs.sbx,'.sbx','.align'),'-mat');

% load frames
% mov = sbxread(dirs.sbx_name,start_frame-1,chunk_size);
nframes = size(RegMat.T,1);
if start_frame+chunk_size > nframes
    mov = PPPack.hf.sbxReadPMT(dirs.sbx,start_frame-1,nframes - start_frame + 1);
else
    mov = PPPack.hf.sbxReadPMT(dirs.sbx,start_frame-1,chunk_size);
end
% squeeze to get rid of other PMTs
mov = squeeze(mov);
% get info variabl
info = PPPack.hf.sbxInfo(dirs.sbx,1);
nframes = info.max_idx+1;
% get correct registration technique
if strcmp(reg_type,'whole pixel')
    ALIGNEDMAT = info.aligned.T;
elseif strcmp(reg_type,'subpixel')
    ALIGNEDMAT = info.aligned.StackReg;
end
if strcmp(reg_type,'whole pixel')
    for i = 1:size(mov,3)
        mov(:,:,i) = circshift(mov(:,:,i),ALIGNEDMAT(start_frame+i-1,:)); % align the image
    end
elseif strcmp(reg_type,'subpixel')
    if start_frame+chunk_size > nframes
        [~,mov] = PPPack.hf.sbxStackRegister(mov,[],[],ALIGNEDMAT(start_frame:end,:),0);
    else
        [~,mov] = PPPack.hf.sbxStackRegister(mov,[],[],ALIGNEDMAT(start_frame:start_frame+chunk_size-1,:),0);
    end
end


