function premasks = remove_double_masks(premasks)
% index back through new masks and make sure that they are not “hollow”
   % i.e, values near centroid have a value > 0 this will prevent using
   % overlapping masks and gives highest priority to early masks with high
   % ica or nmf scores. 
   
ROI_to_chuck = []; 
for curr_ROI1 = 1:length(premasks.masks_keep)-1
    binmask1 = premasks.masks_keep(curr_ROI1).mask > 0;
    mask1 = premasks.masks_keep(curr_ROI1).mask;
    % now iterate through all other masks and subtract the first binary
    % mask from the second binary mask
    for curr_ROI2 = curr_ROI1+1:length(premasks.masks_keep)
        binmask2 = premasks.masks_keep(curr_ROI2).mask > 0;
        mask2 = premasks.masks_keep(curr_ROI2).mask;
        % check if any overlap and if there isnt then continue
        A = binmask1 + binmask2;
        if sum(A(:) > 1) == 0
            continue
        end
        if corr2(mask1,mask2) > 0.4
            B = binmask2-binmask1;
%             stats = regionprops(binmask2);
%             curr_centroid = round(stats.Centroid);
            [~, col] = max(mean(mask2, 1));
            [~, row] = max(mean(mask2, 2));
            if B(row,col) == 0
                ROI_to_chuck = [ROI_to_chuck curr_ROI2];
            end
        end
    end
end
% now remove these from premasks and rebuild masks
% old ROI indices
old_ROIs = find(premasks.ROI_list > 0);
premasks.ROI_list(old_ROIs(ROI_to_chuck)) = 0;
premasks = rmfield(premasks,'masks_keep');
premasks = PPPack.hf.build_masks_for_extraction(premasks);