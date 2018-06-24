function [mov,edges] = sbxRemoveEdges(mov,mov_small)

if nargin < 1
    mov = [];
end

if nargin < 2
    mov_small = [];
end

edges = [70 70 25 25];

% to remove edges
if ~isempty(mov)
    mov = mov(edges(3):end-edges(4)-1,edges(1):end-edges(2)-1,:);
end

% if given small movie get back the edges
if ~isempty(mov_small)
    cols_remove = edges(1);
    rows_remove = edges(3);
    col_buffer = zeros(size(mov_small,1),cols_remove,size(mov_small,3)); 
    % now buffer back in both dimensions
    mov = [col_buffer mov_small col_buffer];
    row_buffer = zeros(rows_remove,size(mov,2),size(mov,3));
    mov = [row_buffer; mov; row_buffer];
end
    

