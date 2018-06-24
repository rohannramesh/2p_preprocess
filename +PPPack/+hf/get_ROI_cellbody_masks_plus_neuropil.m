function signals = get_ROI_cellbody_masks_plus_neuropil(signals)

delete_ROI_tag = 0;

% first lets identify those ROIs that have been labelled non ROIs and
% identify the pixels to not be included in any run
pix_to_not_include_anywhere = zeros(size(signals(1).mask));
ROI_to_eliminate_at_end = [];
ind_cells_remove = 1;
for i = 1:length(signals)
    if isfield(signals(i),'group_number') % to see if that category exists
        if signals(i).group_number == 5 % group 5 is trash
            pix_to_not_include_anywhere = pix_to_not_include_anywhere + ...
                signals(i).mask > 0;
            ROI_to_eliminate_at_end = [ROI_to_eliminate_at_end i];
            signals_remove(ind_cells_remove) = signals(i);
            ind_cells_remove = ind_cells_remove + 1;
        end
    end
end
% now throw the cells to remove from the signals
signals(ROI_to_eliminate_at_end) = [];



%% Make binary masks
for i = 1:length(signals)   
	A = signals(i).mask > 0;
    % don't allow ROI to be in outer most pixel row of FOV
    A(1,:) = 0;
    A(end,:) = 0;
    A(:,1) = 0;
    A(:,end) = 0;
    % don't include any pixels that were eliminate in gui
    A(pix_to_not_include_anywhere > 0) = 0;
    signals(i).binmask = A;    
	signals(i).area = nnz(signals(i).binmask);
end

%% Make overlay cell image

% Dilate cells themselves so have protected ring
all_cells_dilated = zeros(size(signals(1).binmask));
for i = 1:length(signals);
    tmp = signals(i).binmask;
    tmp = imdilate(tmp,[ones(10,10)]); % this is the dilation factor
    all_cells_dilated = all_cells_dilated + tmp;
end
all_cells_dilated = all_cells_dilated > 0;

% Dilation to get neuropil signal
overlap_image = zeros(size(signals(1).binmask));
for i = 1:numel(signals)
    % dilate for building neuropil ring
    inner = imdilate(signals(i).binmask,[ones(15,15)]); % inner edge npil ring
    outer = imdilate(signals(i).binmask,[ones(50,50)]); % outer edge npil ring
    omi = outer - inner; % Have ring around actual axon and have neuropil outside this
    ind_overlap = find(omi == 1 & all_cells_dilated == 1); % Make sure neuropil doesn't overlap with dilated cell
    overlap_image(ind_overlap) = 4;
    omi(ind_overlap) = 0;
    % remove pixels excluded from both cellbodies and neuropil (group 5)
    omi(pix_to_not_include_anywhere > 0) = 0;
    signals(i).neuropil = logical(omi);
    if nnz(signals(i).neuropil)<50
        outer = imdilate(signals(i).binmask,[ones(50,50)]);
        omi = outer - inner; % Have ring around actual axon and have neuropil outside this
        ind_overlap = find(omi == 1 & all_cells_dilated == 1); % Make sure neuropil doesn't overlap with dilated cell
        overlap_image(ind_overlap) = 4;
        omi(ind_overlap) = 0;
        % remove pixels excluded from both cellbodies and neuropil (group 5)
        omi(pix_to_not_include_anywhere > 0) = 0;
        signals(i).neuropil = logical(omi);
    end
end

all_neuropil = zeros(size(signals(1).neuropil)); % Can Neuropil overlap with other neuropil
% figure;
for i = 1:length(signals)
	all_neuropil = all_neuropil + signals(i).neuropil;
end
% imagesc(all_neuropil)

all_cells_overlap = zeros(size(signals(1).binmask));
% To determine overlap image
for i = 1:length(signals)
	all_cells_overlap = all_cells_overlap + signals(i).binmask;
end

% remove all overlapping portions of ROI
ind_overlap = find(all_cells_overlap > 1);
overlap_image = zeros(size(signals(1).binmask));
overlap_image(ind_overlap) = 1;
% No longer dilating image before throwing out
overlap_image_dilated = imdilate(overlap_image,[ones(2,2)]);
overlap_image_dilated = overlap_image;
ind_overlap_dilated = find(overlap_image_dilated > 0);
ind_to_delete = [];
for i = 1:length(signals)
    clear tmp
    tmp = signals(i).binmask;
    tmp(ind_overlap_dilated) = 0;
    signals(i).binmask = tmp;
    % if this removes ROI then remove this from signals
    if max(max(signals(i).binmask)) == 0
        ind_to_delete = [ind_to_delete i];
    end
end
if delete_ROI_tag == 1
    signals(ind_to_delete) = []; 
end

all_cells = zeros(size(signals(1).binmask));
% figure;
for i = 1:length(signals)
	all_cells = all_cells + signals(i).binmask;
end

all_cells = all_cells > 0;
all_cells_undilated = all_cells;



%% Get centroids
% hold on
% 
for i = 1:length(signals)
	[rows, cols] = size(signals(i).binmask);

	y = 1:rows;
	x = 1:cols;

	[X, Y] = meshgrid(x,y);

	signals(i).centroid.x = mean(X(signals(i).binmask==1));
	signals(i).centroid.y = mean(Y(signals(i).binmask==1));
	
	plot(signals(i).centroid.x, signals(i).centroid.y, '.g');
	
end

all_neuropil = all_neuropil > 0;
figure
imagesc(all_neuropil+all_cells_undilated.*2);
for i = 1:numel(signals)
	text(signals(i).centroid.x, signals(i).centroid.y, num2str(i), 'color', [1 0 0],'FontSize',18);
end

%% Adjust ica weights to scale with total pixel number ie will scale as if every pixes had weight of 1

for i = 1:length(signals);
    pix_w = signals(i).mask(signals(i).mask~=0);
    ind_pix = find(signals(i).mask~=0);
    pix_w = pix_w.*(length(pix_w)./sum(pix_w));
    signals(i).mask_norm = zeros(size(signals(i).mask));
    signals(i).mask_norm(ind_pix) = pix_w;
end

%% add an additional mask that is just the neuropil

max_ROIn = length(signals);
allfields = fieldnames(signals(1));
% preallocate all fields
for i = 1:length(allfields)
    signals(max_ROIn+1).(allfields{i}) = [];
end
% Now insert neuropil as mask
signals(max_ROIn+1).mask = all_neuropil > 0;
signals(max_ROIn+1).ica_segment = all_neuropil > 0;
signals(max_ROIn+1).neuropil = zeros(size(all_neuropil));
signals(max_ROIn+1).binmask = all_neuropil > 0;
signals(max_ROIn+1).area = nnz(all_neuropil>0);
% set ica trace = to last ica trace just so don't hit errors - same for
% centroid and 
signals(max_ROIn+1).ica_trace = signals(max_ROIn).ica_trace;
signals(max_ROIn+1).boundary = signals(max_ROIn).boundary;
signals(max_ROIn+1).centroid = signals(max_ROIn).centroid;

