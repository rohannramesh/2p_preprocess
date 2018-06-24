function cellsort = get_ROI_axons_masks_plus_neuropil(cellsort,expt,mouse,run)

delete_ROI_tag = 0
%% Make binary masks
for i = 1:length(cellsort)
%     
	A = cellsort(i).mask > 0;
    % don't allow ROI to be in outer most pixel row of FOV
    A(1,:) = 0;
    A(end,:) = 0;
    A(:,1) = 0;
    A(:,end) = 0;
    cellsort(i).binmask = A;    
	cellsort(i).area = nnz(cellsort(i).binmask);
end

%% Make overlay cell image

% Dilate cells themselves so have protected ring
all_cells_dilated = zeros(size(cellsort(1).binmask));
for i = 1:length(cellsort);
    tmp = cellsort(i).binmask;
    tmp = imdilate(tmp,[ones(4,4)]);
    all_cells_dilated = all_cells_dilated + tmp;
end
all_cells_dilated = all_cells_dilated > 0;

% Dilation to get neuropil signal
overlap_image = zeros(size(cellsort(1).binmask));
for i = 1:numel(cellsort)
    inner = imdilate(cellsort(i).binmask,[ones(4,4)]);
    outer = imdilate(cellsort(i).binmask,[ones(20,20)]);
    omi = outer - inner; % Have ring around actual axon and have neuropil outside this
    ind_overlap = find(omi == 1 & all_cells_dilated == 1); % Make sure neuropil doesn't overlap with dilated cell
    overlap_image(ind_overlap) = 4;
    omi(ind_overlap) = 0;
    cellsort(i).neuropil = logical(omi);
end

all_neuropil = zeros(size(cellsort(1).neuropil)); % Can Neuropil overlap with other neuropil
% figure;
for i = 1:length(cellsort)
	all_neuropil = all_neuropil + cellsort(i).neuropil;
end
% imagesc(all_neuropil)

all_cells_overlap = zeros(size(cellsort(1).binmask));
% To determine overlap image
for i = 1:length(cellsort)
	all_cells_overlap = all_cells_overlap + cellsort(i).binmask;
end

% remove all overlapping portions of ROI
ind_overlap = find(all_cells_overlap > 1);
overlap_image = zeros(size(cellsort(1).binmask));
overlap_image(ind_overlap) = 1;
% No longer dilating image before throwing out
% overlap_image_dilated = imdilate(overlap_image,[ones(3,3)]);
overlap_image_dilated = overlap_image;
ind_overlap_dilated = find(overlap_image_dilated > 0);
ind_to_delete = [];
for i = 1:length(cellsort)
    clear tmp
    tmp = cellsort(i).binmask;
    tmp(ind_overlap_dilated) = 0;
    cellsort(i).binmask = tmp;
    % if this removes ROI then remove this from cellsort
    if max(max(cellsort(i).binmask)) == 0
        ind_to_delete = [ind_to_delete i];
    end
end
if delete_ROI_tag == 1
    cellsort(ind_to_delete) = [];
else
    warning('Not removing ROI that are completely overlapping')
end

all_cells = zeros(size(cellsort(1).binmask));
% figure;
for i = 1:length(cellsort)
	all_cells = all_cells + cellsort(i).binmask;
end

all_cells = all_cells > 0;
all_cells_undilated = all_cells;







%% Get centroids
% hold on
% 
for i = 1:length(cellsort)
	[rows, cols] = size(cellsort(i).binmask);

	y = 1:rows;
	x = 1:cols;

	[X, Y] = meshgrid(x,y);

	cellsort(i).centroid.x = mean(X(cellsort(i).binmask==1));
	cellsort(i).centroid.y = mean(Y(cellsort(i).binmask==1));
	
	plot(cellsort(i).centroid.x, cellsort(i).centroid.y, '.g');
	
end

all_neuropil = all_neuropil > 0;
figure
imagesc(all_neuropil+all_cells_undilated.*2);
for i = 1:numel(cellsort)
	text(cellsort(i).centroid.x, cellsort(i).centroid.y, num2str(i), 'color', [1 1 1],'FontSize',18);
end


%% Adjust ica weights to scale with total pixel number ie will scale as if every pixes had weight of 1

for i = 1:length(cellsort);
    pix_w = cellsort(i).mask(cellsort(i).mask~=0);
    ind_pix = find(cellsort(i).mask~=0);
    pix_w = pix_w.*(length(pix_w)./sum(pix_w));
    cellsort(i).mask_norm = zeros(size(cellsort(i).mask));
    cellsort(i).mask_norm(ind_pix) = pix_w;
end