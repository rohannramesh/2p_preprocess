function premasks = build_masks_for_extraction(premasks)


ROI_to_save = find(premasks.ROI_list > 0);
Group_number = premasks.ROI_list(ROI_to_save);

% Make masks
for i = 1:length(ROI_to_save)
    curr_ROI = ROI_to_save(i);
    % erode if this is a field from the CNN
    if isfield(premasks,'erosions_to_apply')
        curr_mask = premasks.ica(curr_ROI).filter;
        curr_mask = PPPack.hf.get_erosion_mask(curr_mask,premasks.erosions_from_CNN(curr_ROI));
    else
        curr_mask = premasks.ica(curr_ROI).filter;
    end
    
    filtered = PPPack.hf.imgGaussBlur(curr_mask, 1);
    [contourCoords,~] = contour(filtered, [1,1]*(mean(filtered(:))+premasks.countour_sd*std(filtered(:))));

    % Clean up contour:
    keep = abs(diff(contourCoords(1, :))) < 50;
    keep = keep .* (abs(diff(contourCoords(2, :))) < 50);
    keep = [keep 1]; % Compensate for element lost by diff();	
    keep = logical(keep);
    cleanX = contourCoords(1, keep);
    cleanY = contourCoords(2, keep);
    premasks.ica(curr_ROI).contour = [cleanX; cleanY];		
end
								

for i = 1:length(ROI_to_save)
    curr_ROI = ROI_to_save(i);
        if isfield(premasks,'erosions_to_apply')
            curr_mask = premasks.ica(curr_ROI).filter;
            curr_mask = PPPack.hf.get_erosion_mask(curr_mask,premasks.erosions_from_CNN(curr_ROI));
        else
            curr_mask = premasks.ica(curr_ROI).filter;
        end
		currentICA.mask = curr_mask;
		currentICA.boundary = premasks.ica(curr_ROI).contour;
		currentICA.trace = [];
        currentICA.ica_trace = premasks.ica(curr_ROI).trace; % added RR
        currentICA.group_number = Group_number(i);
		premasks.masks_keep(i) = currentICA;
end