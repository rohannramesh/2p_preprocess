function curr_norm_hist = get_hist_of_tc(premasks,curr_ROI,bins_use)
    curr_trace = premasks.ica(curr_ROI).trace;
    curr_trace = curr_trace - median(curr_trace);
    curr_trace = curr_trace./max(curr_trace);
    % binsizes
   curr_norm_hist = hist(curr_trace,bins_use);
    curr_norm_hist = curr_norm_hist./sum(curr_norm_hist);