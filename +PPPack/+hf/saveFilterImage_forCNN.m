function [cx, cy, ran] = saveFilterImage_forCNN(path, im, minim)
% save image of each filter for CNN testing

    if nargin < 3, minim = 20; end
    
    A = im(:);
    B = PPPack.hf.extract_patch(A,size(im),[minim minim]);
    smallim = B;
    if size(smallim,1) ~= minim || size(smallim,2) ~= minim
        return
    end
    imwrite(smallim, path);
    
    


end

