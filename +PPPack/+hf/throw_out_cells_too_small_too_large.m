function premasks = throw_out_cells_too_small_too_large(premasks,min_ROI_size,max_ROI_size);
% throw out ROIs too small or too large

ROIs_to_chuck = [];
for i = 1:length(premasks.ica)
    curr_mask = premasks.ica(i).filter > 0;
    if sum(curr_mask(:)) <= min_ROI_size | sum(curr_mask(:)) >= max_ROI_size
        ROIs_to_chuck = [ROIs_to_chuck; i];
    end
end

% now update variable
premasks.suggested_ROI_keep(ROIs_to_chuck) = [];
premasks.suggested_ROI_nonkeep(ROIs_to_chuck) = [];
premasks.ROIs_selected_from_orig(ROIs_to_chuck) = [];
premasks.ROI_list(ROIs_to_chuck) = [];
premasks.ica(ROIs_to_chuck) = [];
premasks.Sorting_metric(ROIs_to_chuck) = [];
premasks.Weights.Spatial(ROIs_to_chuck) = [];
premasks.Weights.Temporal(ROIs_to_chuck) = [];


