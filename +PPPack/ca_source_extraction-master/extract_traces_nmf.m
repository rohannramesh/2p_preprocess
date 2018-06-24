function traces = extract_traces_nmf(Masks,data,chunk_size)

if nargin < 3
    chunk_size = 1e4;
end

sizY = size(data);
T = sizY(end);

[d1a,d2a] = size(Masks);

d1y = prod(sizY(1:end-1));
d2y = T;

traces = zeros(d2a,d2y);
chunks = 1:chunk_size:d1a; corenum = 10;
nchunks = numel(chunks); % number of indexes
chunk_idx = arrayfun(@(i) chunks(i:min((i+corenum-1), nchunks)), ...
    1:corenum:nchunks, 'UniformOutput', false); % indices of the patches in each batch
for i = 1:numel(chunk_idx)
    batch2run = chunk_idx{i};
    % predetermine chunks to pass to parfor
    clear tmp_mov_chunks
    for ii = 1:numel(batch2run)
        idx = batch2run(ii);
        tmp_mov_chunks{ii} = data(idx:min(idx+chunk_size-1,d1a),:);
    end
    parfor ii = 1:numel(batch2run)
%     for ii = 1:numel(batch2run)
        idx = batch2run(ii);
%         tempY = data(idx:min(idx+chunk_size-1,d1a),:);
        tempY = tmp_mov_chunks{ii};
        tempA = Masks(idx:min(idx+chunk_size-1,d1a),:);
        AYt{ii} = tempA'*double(tempY);
    end
    traces = traces + sum(cat(3, AYt{:}),3);
    clear AYt
    if mod(i, 20) == 0
        fprintf('%2.1f%% of chunks completed \n', i*100/numel(chunk_idx));
    end
end