function ROIvars = add_tc_to_ROIvars(data,A,b,ROIvars);

sizY = size(data,'Y');
A = double(A);
% get tc
if ~isempty(b)
    AY = mm_fun([A,double(b)],data);
    K = size(A,2);
    bY = AY(K+1:end,:);
    AY = AY(1:K,:);
    % now normalize by pixel number
    tc = bsxfun(@rdivide,AY,sum(A,1)'); % 4 is bc of spatial downsampling
    tc_b = bY./(sum(b)');
    % now npil correct
    tc = bsxfun(@minus,tc,tc_b);
else
    AY = mm_fun([A],data);
    K = size(A,2);
%     bY = AY(K+1:end,:);
    AY = AY(1:K,:);
    % now normalize by pixel number
    tc = bsxfun(@rdivide,AY,sum(A,1)'); % 4 is bc of spatial downsampling
end
ROIvars.tc = tc;


