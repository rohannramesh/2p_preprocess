function signals = deconvolveRemovingNANRows(signals)
%DECONVOLVEREMOVINGNANROWS Is required in the rare case where one ROI is
%   not found across an entire run. It removes NAN rows from the dff
%   signal, submits the rest for deconvolution, and correctly pads upon
%   return.


    % make dff variable
    nROIs = length(signals.timecourse);
    dff = nan(nROIs,length(signals.timecourse(1).dff_sliding));
    for curr_ROI = 1:nROIs
        dff(curr_ROI,:) = signals.timecourse(curr_ROI).dff_sliding;
%         dff(curr_ROI,:) = signals.timecourse(curr_ROI).dff_axon;
    end

    nanrows = zeros(1, size(dff, 1));
    for c = 1:size(dff, 1)
        if sum(isnan(dff(c, :))) > size(dff, 2)/2
            nanrows(c) = 1;
        end
    end

    if sum(nanrows) == 0
        decon = PPPack.hf.calc_constrained_foopsi(dff);
        deconvolved = decon.deconvInSpikes;
    else
        submission = dff(nanrows == 0, :);
        
        decon = PPPack.hf.calc_constrained_foopsi(submission);
        out = decon.deconvInSpikes;
        
        deconvolved = NaN(size(dff, 1), size(dff, 2));
        deconvolved(nanrows == 0, :) = out;
    end
    
    % put into signals file
    for curr_ROI = 1:nROIs
        signals.timecourse(curr_ROI).deconvolved = deconvolved(curr_ROI,:);
    end
    
end

