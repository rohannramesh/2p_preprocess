function premasks = run_PCAICA(movpath,temp_dir,mov)

if nargin < 2
	temp_dir = [];
end

if nargin < 3;
    mov = [];
end


nPCs = 1000;
flims = [];
dsamp = [];
segment = true; fprintf('Segmenting: %s\n', segment);
mu = 0.5; fprintf('Mu (0 = spatial): %d\n', mu); 


[mixedsig, mixedfilters, CovEvals, covtrace, meanimg, ~] = ...
    PPPack.PCAICA.CellsortPCA_RR(mov, flims, nPCs, dsamp, temp_dir, []);






% 2a. Choose PCs

% [PCuse] = CellsortChoosePCs(movpath, mixedfilters);
 
% 2b. Plot PC spectrum
% figure;
last_pc_above_noise = PPPack.PCAICA.CellsortPlotPCspectrum_RR(mov, CovEvals, 1:nPCs, 1);
% last_pc_above_noise = nPCs;
% Discard first three PCs because they are likely to contain full-field
% effects:
PCuse = 4:last_pc_above_noise; 


% 3a. ICA
nIC = length(PCuse);

[ica_sig, ica_filters, ica_A, numiter] = ...
	PPPack.PCAICA.CellsortICA(mixedsig, mixedfilters, CovEvals, PCuse, mu, nIC);

% Normalise ICA_filters such that they are not negative...may not make
% super much sense but I got this from the original CellsortSegmentation
% function:
ica_filters = (ica_filters - mean(ica_filters(:)))/abs(std(ica_filters(:)));


%% 4. segment
if segment

	smwidth		= 2;	% Standard deviation of Gaussian smoothing kernel (pixels)
	thresh		= 2;	% Threshold for spatial filters (standard deviations) - used to be 5 - RR
	arealims	= 25;	% 1 (or 2)-element vector specifying the minimum (and maximum) area
	plotting	= 0;
	[ica_segments, segmentlabel, segcentroid] = ...
		PPPack.PCAICA.CellsortSegmentation(ica_filters, smwidth, thresh, arealims, plotting);

	% Get segment time series
	segment_sig = PPPack.PCAICA.CellsortApplyFilter_nopath(mov, ica_segments, flims, meanimg, 1);
    
else
	ica_segments = ica_filters;
	segment_sig = ica_sig;
	segcentroid = zeros(size(ica_segments, 1), 1);
end

%%  make variables global and call display GUI
if nargout < 1
	global premasks
end

% Initialize premasks variable and put masks inside
premasks = [];
premasks.countour_sd = 5; % Set the stringency with which ROIs are generated from ICs
for i = 1:size(ica_segments, 1)
    premasks.ica(i).filter = squeeze(ica_segments(i,:,:));
    premasks.ica(i).trace = segment_sig(i,:);
	premasks.ica(i).centroid = segcentroid(i);
    premasks.ica(i).roiNum = []; %Initialize roinum field.
end


% order the premasks output by the signal to noise ratio from SugDog
try premasks_backup = premasks;
    premasks2 = premasks_backup;
    for i = 1:length(premasks2.ica)
        curr_trace = premasks2.ica(i).trace;
        signal = PPPack.hf.dffsnr(curr_trace);
        signal_per_cell(i) = signal;
    end
    [A,sort_ind] = sort(signal_per_cell,'descend');
    for i = 1:length(premasks2.ica)
        premasks2.ica(i).filter = premasks_backup.ica(sort_ind(i)).filter;
        premasks2.ica(i).trace = premasks_backup.ica(sort_ind(i)).trace;
        premasks2.ica(i).centroid = premasks_backup.ica(sort_ind(i)).centroid;
        premasks2.ica(i).roiNum = premasks_backup.ica(sort_ind(i)).roiNum;
    end
    premasks = premasks2;
catch err
end

premasks.movm = meanimg; %Mean image
premasks.fName = movpath;

% get pixel by pixel corr image
try mov_tds = PPPack.hf.group_z_project(mov,100);
    tic
    A = PPPack.hf.CrossCorrImage(mov_tds);
    premasks.movcorr = A;
    toc
catch err
end

% additionally set up other variables of interest
for i = 1:length(premasks.ica)
    [premasks.ica(i).idx(:,1) premasks.ica(i).idx(:,2)] = find(single(premasks.ica(i).filter>0));
end

% Mask variable
premasks.masks = zeros(size(premasks.movm));
premasks.masks_with_group = zeros(size(premasks.movm));

% Set up ROI_list
premasks.ROI_list = zeros(1,length(premasks.ica));

