function icaguidata = run_PCAICA_Axon(movpath,temp_dir,mov)

if nargin < 2
	temp_dir = [];
end


nPCs = 1000;
flims = [];
dsamp = [];
segment = true; fprintf('Segmenting: %s\n', segment);
mu = 0.9; fprintf('Mu (0 = spatial): %d\n', mu); % Default was 0.5... 


[mixedsig, mixedfilters, CovEvals, covtrace, meanimg, ~] = ...
    PPPack.PCAICA.CellsortPCA_RR(mov, flims, nPCs, dsamp, temp_dir, [], 1);




% 2a. Choose PCs
last_pc_above_noise = PPPack.PCAICA.CellsortPlotPCspectrum_RR(mov, CovEvals, 1:nPCs, 1);
% Discard first three PCs because they are likely to contain full-field
% effects:
PCuse = 4:last_pc_above_noise; % used to be 4:



% 3a. ICA
nIC = length(PCuse);

[ica_sig, ica_filters, ica_A, numiter] = ...
	PPPack.PCAICA.CellsortICA(mixedsig, mixedfilters, CovEvals, PCuse, mu, nIC);

% Normalise ICA_filters such that they are not negative...may not make
% super much sense but I got this from the original CellsortSegmentation
% function:
% ica_filters = (ica_filters - mean(ica_filters(:)))/abs(std(ica_filters(:)));


%% 4. segment
if segment

	smwidth		= 1;	% Standard deviation of Gaussian smoothing kernel (pixels)
	thresh		= 2;	% Threshold for spatial filters (standard deviations) - used to be 5 - RR
	arealims	= 50;	% 1 (or 2)-element vector specifying the minimum (and maximum) area
	plotting	= 0;
	[ica_segments, segmentlabel, segcentroid] = ...
		PPPack.PCAICA.CellsortSegmentation_axon(ica_filters, smwidth, thresh, arealims, plotting);

	% Get segment time series
	segment_sig = PPPack.PCAICA.CellsortApplyFilter_nopath(mov, ica_segments, flims, meanimg, 1);

    
else
	ica_segments = ica_filters;
	segment_sig = ica_sig;
	segcentroid = zeros(size(ica_segments, 1), 1);
end

%% Create relevant variables
if nargout < 1
	global icaguidata
end

% Initialize icaguidata variable and put masks inside
icaguidata = [];
icaguidata.countour_sd = 5; % Set the stringency with which ROIs are generated from ICs
for i = 1:size(ica_segments, 1)
    icaguidata.ica(i).filter = squeeze(ica_segments(i,:,:));
    icaguidata.ica(i).trace = segment_sig(i,:);
	icaguidata.ica(i).centroid = segcentroid(i);
    icaguidata.ica(i).roiNum = []; %Initialize roinum field.
end


% order the icaguidata output by the signal to noise ratio from SugDog
try icaguidata_backup = icaguidata;
    icaguidata2 = icaguidata_backup;
    for i = 1:length(icaguidata2.ica)
        curr_trace = icaguidata2.ica(i).trace;
        signal = PPPack.hf.dffsnr(curr_trace);
        signal_per_cell(i) = signal;
    end
    [A,sort_ind] = sort(signal_per_cell,'descend');
    for i = 1:length(icaguidata2.ica)
        icaguidata2.ica(i).filter = icaguidata_backup.ica(sort_ind(i)).filter;
        icaguidata2.ica(i).trace = icaguidata_backup.ica(sort_ind(i)).trace;
        icaguidata2.ica(i).centroid = icaguidata_backup.ica(sort_ind(i)).centroid;
        icaguidata2.ica(i).roiNum = icaguidata_backup.ica(sort_ind(i)).roiNum;
    end
    icaguidata = icaguidata2;
catch err
end

icaguidata.movm = meanimg; %Mean image
icaguidata.fName = movpath;

% get pixel by pixel corr image
try mov_tds = PPPack.hf.group_z_project(mov,100);
    tic
    A = PPPack.hf.CrossCorrImage(mov_tds);
    icaguidata.movcorr = A;
    toc
catch err
end

% additionally set up other variables of interest
for i = 1:length(icaguidata.ica)
    tmpR = find(single(icaguidata.ica(i).filter>0));
    [icaguidata.ica(i).idx(:,1) icaguidata.ica(i).idx(:,2)] = find(single(icaguidata.ica(i).filter>0));
    if isempty(tmpR)
        [icaguidata.ica(i).idx(:,1)] = [1];
        [icaguidata.ica(i).idx(:,2)] = [1];
    end
end

% Mask variable
icaguidata.masks = zeros(size(icaguidata.movm));
icaguidata.masks_with_group = zeros(size(icaguidata.movm));

% Set up ROI_list
icaguidata.ROI_list = zeros(1,length(icaguidata.ica));



