function [new_vec] = downsample_vector(vector,ds_val,ds_approach)

if nargin < 3
    ds_approach = 'mean';
end

% make sure right orientation
if size(vector,2) == 1
    vector = vector';
end

sz = size(vector);
% groupZ_value = sz(3)/slice_number;
slice_number = sz(2)/ds_val;
vector2 = reshape(vector,sz(1),ds_val,slice_number);
if strcmp(ds_approach,'sum')
    new_vec = sum(vector2,2);
elseif strcmp(ds_approach,'mean')
    new_vec = mean(vector2,2);
elseif strcmp(ds_approach,'max')
    new_vec = max(vector2,[],2);
end
std_zproj = std(vector2,[],3);
SE_zproj = (std(vector2,[],3))./slice_number;
new_vec = squeeze(new_vec); % Group Z project by appropriate number
std_zproj = squeeze(std_zproj);
SE_zproj = squeeze(SE_zproj);
