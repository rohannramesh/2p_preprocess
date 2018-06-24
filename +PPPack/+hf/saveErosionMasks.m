function [maskinds, erosions] = saveErosionMasks(path, filter, erosions, connected)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 4, connected = true; end
    
    % Sort erosions so that we can build up masks
    erosions = sort(erosions);
    erosions = [erosions 2.0];
    
    % Get only those values involved in mask
%     flatnan = filter(:);
%     vals = sort(flatnan(flatnan ~= 0));
%     
%     thresh = zeros(1, length(erosions));
    masks = cell(1, length(erosions));
    maskinds = cell(1, length(erosions));
    
    % Write the beginning of the json file, if desired
    if ~isempty(path)
        fp = fopen(path, 'w');
        fprintf(fp, '{"erosions":[');
        for i = 1:length(erosions)
            if i > 1, fprintf(fp, ','); end
            fprintf(fp, '%.2f', erosions(i));
        end
        fprintf(fp, '],\n"masks":[');
    end
    
%     mn = abs(min(vals)) + 1;
%     filter(filter == 0) = NaN;
%     filter = filter + mn;
%     filter(isnan(filter)) = 0;
    
    for i = 1:length(erosions)
        % Find the threshold and binarize the mask
%         if i == length(erosions)
%             thresh(i) = 0.0000001;
%         else
%             n = round(erosions(i)*length(vals));
%             if n > length(vals), n = length(vals); end
%             if n < 1, n = 1; end
%             thresh(i) = vals(end - n + 1);
%         end
%             
%         mask = filter(:, :);
%         mask(mask < thresh(i)) = 0;
%         mask(mask > 0) = 1;
%         
%         % Only pay attention to the smallest ROI
%         if connected && i < length(erosions)
%             cc = bwconncomp(mask);
%             if cc.NumObjects > 1
%                 [~, mxpos] = max(filter(:));
%                 [y, x] = ind2sub(size(filter), mxpos);
%                 maskfill = imfill(logical(mask - 1), [y x]);
%                 mask = logical(mask) & maskfill;
%             end
%         end
%         
%         % Fill in 'hole' pixels
%         mask = imfill(logical(mask), 'holes');

        mask = getErosionMask(filter, erosions(i), connected);
        
        for j = 1:i-1
            mask = mask - masks{j};
        end
        masks{i} = mask;
        
        % Make a list of x, y indices that make up the mask
        for x = 1:size(mask, 2)
            for y = 1:size(mask, 1)
                if masks{i}(y, x) > 0
                    maskinds{i} = [maskinds{i} x y];
                end
            end
        end
        
        % Write the indices, if desired
        if ~isempty(path)
            if i > 1, fprintf(fp, ','); end
            fprintf(fp, '[');
            for j = 1:length(maskinds{i})
                if j > 1, fprintf(fp, ','); end
                fprintf(fp, '%i', maskinds{i}(j));
            end
            fprintf(fp, ']');
        end
    end
    
    % Write the end and close the file, if desired
    if ~isempty(path)
        fprintf(fp, ']}'); 
        fclose(fp);
    end
end

