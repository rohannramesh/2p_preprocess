function plot_overlay_image_masks(curr_dir,obj,signals)


start_frame = 1000;
frames_to_load = 1000; 

% first lets load second 1000 reg frames
mov = PPPack.hf.sbxLoadReg(curr_dir,obj.PreProcessingParameters.registration,start_frame,frames_to_load);

% make masks image
all_masks = zeros(size(mov,1),size(mov,2));
nROIs = length(signals.timecourse)-1;
for curr_ROI = 1:nROIs
    ind_ROI = find(signals.maskNumber == curr_ROI);
    all_masks(ind_ROI) = 1;
end
    
%%
% now lets overlay mean image
figure('units','normalized','outerposition',[0 0 1 1])
img = nanmean(mov,3);
img(img<0) = 0;
img(isnan(img)) = 0;
img = sqrt(img);
img = img/max(img(:));
img = adapthisteq(img);
tmpRRR = median(double(nonzeros(img)));
imagesc(img,[0 tmpRRR*2.5]);
colormap('gray')
axis image
axis off


mask_group_colors = [0 255 255;... % light blue
                        255 140 0;... % orange
                        50 205 50;... % lime green
                        128 0 0;... % dark salmon
                        0 0 128;... % add many lines of navy, so that if there are too many overlapping regions the GUI doesnt crash
                        0 0 128;... % again, navy
                        0 0 128;... % navy
                        0 0 128;... % naby
                        0 0 128]./255;  % navy
                    
effective_group_number = unique(all_masks(:));
effective_group_number(effective_group_number == 0) = [];

% plot masks that have been accepted
if ~isempty(effective_group_number)
    for curr_group = effective_group_number'
    %     mask_img = icaguidata.masks;
    %     mask_img = mask_img > 0;
        mask_img = all_masks == curr_group;
    %     mask_color = cat(3, 0*ones(size(mask_img)), 255*ones(size(mask_img)), 255*ones(size(mask_img)))./255;
        mask_color = cat(3, mask_group_colors(curr_group,1)*ones(size(mask_img)), ...
            mask_group_colors(curr_group,2)*ones(size(mask_img)), ...
            mask_group_colors(curr_group,3)*ones(size(mask_img)));
        hold on
        m = imagesc(mask_color);
        hold off
        transparency_factor = 0.1;
        set(m,'AlphaData',mask_img.*transparency_factor)
    end
end

% add on numbers
hold on
for curr_ROI = 1:nROIs
    curr_mask = signals.maskNumber == curr_ROI;
    stats = regionprops(curr_mask);
    try text(stats(1).Centroid(1)+3,stats(1).Centroid(2)+3, num2str(curr_ROI),...
        'Color', [0 0 0],...
        'FontSize', 8,...
        'FontWeight', 'bold',...
        'HorizontalAlignment', 'center');
    catch err
    end
end
% truesize()
if strcmp(obj.PreProcessingParameters.ROI_algorithm,'NMF')
    PPPack.hf.screen2tiff([curr_dir.path '\FOV+ROI_nmf'])
else
    PPPack.hf.screen2tiff([curr_dir.path '\FOV+ROI'])
end
