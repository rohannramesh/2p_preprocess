function PreProcessingParameters = get_PreProcessingParameters(PreProcessInputs)
% these are the default settings for each user that they can set
% individually or use commonly set default settings

if strcmp(username,'rramesh1')
    defaults.rig                          = 'Hutch';
    defaults.registration                 = 'subpixel'; % other option is 'whole pixel'
    defaults.write_first_last_500         = 1; % do you want to only write the first and last 500 frames to check registration - make 0 to write entire mov as tiff stack
    defaults.s_before                     = 1; % s before stim onset for checking visually evoked images
    defaults.s_after                      = 2; % s after stim onset for checking visually evoked images
    defaults.clean_data                   = 'no';
    defaults.ROI_selection_process        = 'GUI'; % options are GUI or CNN
    defaults.ROI_type                     = 'Cell'; % options are Cell or Axon
    defaults.npil_correction              = 'Normal'; % other option is Weighted
    defaults.ROI_algorithm                = 'PCA/ICA'; % other option is NMF
    defaults.Cell_Width                   = 2.5; % default half width
    defaults.Patch_Size                   = [30, 30]; % default patch size
    defaults.start_time                   = -1; % time to start analysis - -1 means immediately
    
    
elseif strcmp(username,'cburgess')
    defaults.rig                          = 'Hutch';
    defaults.registration                 = 'subpixel'; % other option is 'whole pixel'
    defaults.write_first_last_500         = 1; % do you want to only write the first and last 500 frames to check registration
    defaults.s_before                     = 1; % s before stim onset for checking visually evoked images
    defaults.s_after                      = 2; % s after stim onset for checking visually evoked images
    defaults.clean_data                   = 'no';
    defaults.ROI_selection_process        = 'GUI'; % options are GUI or CNN
    defaults.ROI_type                     = 'Cell'; % options are Cell or Axon
    defaults.npil_correction              = 'Normal'; % other option is Weighted
    defaults.ROI_algorithm                = 'PCA/ICA'; % other option is NMF
    defaults.Cell_Width                   = 2.5; % default half width
    defaults.Patch_Size                   = [30, 30]; % default patch size    
    defaults.start_time                   = -1; % time to start analysis - -1 means immediately
    
elseif strcmp(username,'kmcguir2') || strcmp(username,'kmcguire')
    defaults.rig                          = 'Hutch';
    defaults.registration                 = 'subpixel'; % other option is 'whole pixel'
    defaults.write_first_last_500         = 1; % do you want to only write the first and last 500 frames to check registration
    defaults.s_before                     = 1; % s before stim onset for checking visually evoked images
    defaults.s_after                      = 2; % s after stim onset for checking visually evoked images
    defaults.clean_data                   = 'no';
    defaults.ROI_selection_process        = 'CNN'; % options are GUI or CNN
    defaults.ROI_type                     = 'Cell'; % options are Cell or Axon
    defaults.npil_correction              = 'Normal'; % other option is Weighted
    defaults.ROI_algorithm                = 'PCA/ICA'; % other option is NMF
    defaults.Cell_Width                   = 2.5; % default half width
    defaults.Patch_Size                   = [30, 30]; % default patch size    
    defaults.start_time                   = -1; % time to start analysis - -1 means immediately
    
else %%%%%%% FOR GENERAL USERS
    defaults.rig                          = 'Starsky';
    defaults.registration                 = 'subpixel'; % other option is 'whole pixel'
    defaults.write_first_last_500         = 1; % do you want to only write the first and last 500 frames to check registration
    defaults.s_before                     = 1; % s before stim onset for checking visually evoked images
    defaults.s_after                      = 2; % s after stim onset for checking visually evoked images
    defaults.clean_data                   = 'no';
    defaults.ROI_selection_process        = 'GUI'; % options are GUI or CNN
    defaults.ROI_type                     = 'Cell'; % options are Cell or Axon
    defaults.npil_correction              = 'Normal'; % other option is Weighted
    defaults.ROI_algorithm                = 'PCA/ICA'; % other option is NMF
    defaults.Cell_Width                   = 2.5; % default half width
    defaults.Patch_Size                   = [30, 30]; % default patch size    
    defaults.start_time                   = -1; % time to start analysis - -1 means immediately
end

% first lets set our variable equivalent to the default and then replace
% with any specifics entered by the user
PreProcessingParameters = defaults;
if ~isempty(PreProcessInputs)
    curr_fieldnames_choice = fieldnames(PreProcessInputs);
    curr_fieldnames_default = fieldnames(defaults);
    for i = 1:length(curr_fieldnames_choice)
        if sum(ismember(curr_fieldnames_default,curr_fieldnames_choice{i})) == 0
            error(['PreProcessing parameter input not found: ' curr_fieldnames_choice{i}]);
        end
        PreProcessingParameters.(curr_fieldnames_choice{i}) = PreProcessInputs.(curr_fieldnames_choice{i});
    end
end