%% Information about running the Pre-processing pipeline

% This script will run the entire Pre-Processing package for your 2-photon
% imaging data set. Here, I will describe briefly the approach taken at
% each step, the parameters that should be set depending on your particular
% imaging setup, and decisions that need to be made in terms of how you
% want to preprocess your data. Make sure you read the GitHub ReadME or 
% the Setup section below for proper initialization.
% Within the PPPack (Pre-Process Package)
% there are three main folders: the PreProcess2P class, the helper
% functions package (hf), the constrained non-negative matrix
% factorization (NMF) package, and the PCA/ICA package. These last two folders
% contain the functions and pre-trained convolutional neural networks needed to
% identify putative ROIs and for ROI selection

%% Setup

% All code has been tested and optimized in Matlab 2015b+. It is unclear if
% all functions will work the same in earlier versions. 
% Dependencies: Statistics and Machine Learning Toolbox, Image Processing
% toolbox
% Recommended but not necessary: Parallel compution toolbox (for large data
% sets) and Signal processing toolbox

% Install python v 3.6 and NumPy, Keras, and Tensorflow if you would
% like to use the pre-trained convolutional neural networks for ROI selection and make sure
% Matlab acknowledges Python's existence using pyversion.

% From within +PPPack move ca_source_extraction_master and glmnet_matlab into your path.
% Also add both ROI selection GUI's for both PCA/ICA and NMF to
% your path (both .fig and m file).
% All folders with a + in front of it are protected and can not be added to
% the path. 
% For deconvolution, we are using an algorithm that requires the CVX
% library, download here:  http://cvxr.com/cvx/download/

%% Preprocess pipeling

% PreProcessing Steps (for details see each cell below)
% 1. Define your PreProcessing object
% 2. Save object
% 3. Register data
% 4. Build frame metadata and trial information - (optional) if using MonkeyLogic behavioral training setup
% (Asaad and Esakandar, 2008)
% 5. ROI identification - PCA/ICA or CNMF
% 6. ROI extraction - pre-trained convolutional neural network or ROI
% selection GUI (for semi-automated manual selection)
% 7. Signal fluorescence trace extraction (neuropil correction and dFF
% calculation)
% 8. Run Generalized Linear Model (optional)

%% General imaging protocol used

% Each 2P imaging session can contain multiple imaging runs. In the Andermann
% lab we collect ~4 30 min imaging runs at a single plane. A separate
% object should be created for each new plane imaged. Within a given
% imaging plane there will be multiple active neurons. The goal of this
% preprocessing pipeline is to extract the activity of these neurons for
% all imaging runs on a given day. Our lab uses the relatively common
% Scanbox imaging setup and so within the folder for an imaging run there
% will be the following files:
%     1. .sbx - raw data from that imaging run (still have to make flexible
%               for tiff setup)
%     2. .ephys - nidaq data if aligning animal behavior to 2P frames
%     3. .mat - info file for 2P imaging run (paired with .sbx file)
%     4. _eye.mat - file with eye tracking data
%     5. .quadrature - information about animal's running
% 
% the pre-processing files will build the following files (load each
% filetype as a common .mat file- format used is for easy file identification)
%     1. .align - file containing the subpixel shifts for each plane to 
%                 register each run to a common target
%     2. .f2p and .TrialVar - file for information about stimulus and behavior on each frame
%     3. .ica or .nmf - output from PCA/ICA or CNMF ROI identification algorithm respecttively
%     4. .icamasks or .nmfmasks - output after ROI selection process (GUI or CNN)
%     5. .npilmasks - masks for each neuropil mask for each ROI identified
%     6. .signalsica or .signalsnmf - signal variable containing the extracted timecourses
% 
% Folder organization we used:
% All of the data for each run was saved in a single folder, with multiple
% folders for multiple runs for each day.
% An example path to a single run folder:
% 'S:\twophoton_data\2photon\scan\OA27\170524_OA27\170524_OA27_run4'
% The important considerations are Base\MouseName\Date\Runs_to_analyze
% We define base in PPPack.hf.sbxScanBase and everything else in PPPack.hf.sbxDir
% From personal experience we recommend this structure, but feel free to
% adjust accordingly.



%% Define PreProcessing Object

% the object (myObj) can be initialized one of two ways. First you can pass
% empty square brackets (i.e. PPPack.PreProcess2P([])). This will pop up a
% dialogue box that allows you to enter the MouseName, the Dates, and the
% Runs to analyze. The other option allows you to pass the MouseName,
% Dates, and Runs_to_analyze directly
% (PPPack.PreProcess2P('RR1','170524',2:4)).
% PreProcessing parameters to choose during creation of object:
% defaults for each user can be set in PPPack.hf.get_PreProcessingParameters

%     defaults.rig                          = 'Hutch'; % Rig name if multiple 2P rigs exist
%     defaults.registration                 = 'subpixel'; % other option is 'whole pixel'
%     defaults.write_first_last_500         = 1; % do you want to only write the first and last 500 frames to check registration - make 0 to write entire mov as tiff stack
%     defaults.ROI_selection_process        = 'GUI'; % options are GUI (manual) or CNN (automated)
%     defaults.ROI_type                     = 'Cell'; % options are Cell or Axon
%     defaults.npil_correction              = 'Normal'; % other option is Weighted
%     defaults.ROI_algorithm                = 'PCA/ICA'; % other option is NMF
%     defaults.Cell_Width                   = 2.5; % default half width - only used for CNMF
%     defaults.Patch_Size                   = [30, 30]; % default patch size
%     defaults.s_before                     = 1; % s before stim onset for checking visually evoked images
%     defaults.s_after                      = 2; % s after stim onset for checking visually evoked images
%     defaults.clean_data                   = 'no';

% Alternatively you can create a variable such as described below:
    % PreProcessOptions.registration = 'subpixel';
    % PreProcessOptions.clean_data   = 'no';
    % PreProcessOptions.npil_correction     = 'weighted';
    % PreProcessOptions.ROI_algorithm = 'NMF';
    % PreProcessOptions.rig = 'Starsky';
% The other options will be autopopulated based on your personal settings
% or the default settings and can be passed in to create the object as
% follows: PPPack.PreProcess2P('RR1','170524',2:4,PreProcessOptions)

myObj = PPPack.PreProcess2P('MouseName','Date',Runs_to_use)


%% Save object.
% Object will be saved in the folder for a given day (depends on data
% structure as defined by PPPack.hf.sbxScanBase

myObj.save2PPPObj;

%% Register data
% Each frame on each run will be registered to a target image from the
% first run included in the object using fast, efficient sub-pixel
% approaches (Bonin et al., 2011). This will save a .align file within the
% folder for each run containing the amount to shift each frame for proper
% alignment

myObj.sbxRegister;

%% Test registration via tiff stack visualization
% Write a tiff stack that allows you to quickly check if the registration
% has worked as expected. Per run this will subselect the first and last
% 500 frames for visualization

myObj.sbxTestRegistrationImageStack;

%% Build metadata and trial variable files (if using MonkeyLogic as training rig )
% These two functions will align all visual stimulus and behavior data to
% the appropriate 2p framerate and create metadata variables so you have
% information about exactly what the animal is seeing and doing for every
% single data frame

myObj.sbxframe2pmetadata;
myObj.sbxTrialInfo;

%% Run ROI identification algorithm
% Either identify putative regions of interest (ROIs) using either the
% PCA/ICA algorithm from Mukamel and Schnitzer (2009) or the CNMF algorithm
% (Pnevmatikakis et al., 2016). This step will build either a .ica or .nmf
% file for future ROI selection. These files contain a variable called
% premasks which contains all possible ROIs identified by the respective
% algorithm, their weights, and a timecourse associated with them. To
% generate these masks we spatially downsample by 2 and temporally
% downsample by 5 for time-saving purposes

if strcmp(myObj.PreProcessingParameters.ROI_algorithm,'PCA/ICA')
    myObj.PCAICA_ROI_identification;
elseif strcmp(myObj.PreProcessingParameters.ROI_algorithm,'NMF')
    myObj.CNMF_ROI_identification;
end
    
%% ROI selection
% Choose which ROIs to keep from either ROI identification algorithm using
% either a pre-trained convolutional neural network, or manual selection
% using a custom-designed Matlab GUI

% Useful GUI hotkeys
% n - next, b - back, k - keep

% The output of this step will be a variable called .icamasks or .nmfmasks
% which has a subtructure within premasks containing the manual or CNN selected ROIs to
% keep. 

if strcmp(myObj.PreProcessingParameters.ROI_selection_process,'GUI')
    myObj.run_ROI_selection_gui;
elseif strcmp(myObj.PreProcessingParameters.ROI_selection_process,'CNN')  
    if strcmp(myObj.PreProcessingParameters.ROI_algorithm,'PCA/ICA')
        myObj.ROI_selection_CNN_PCAICA()
    elseif strcmp(myObj.PreProcessingParameters.ROI_algorithm,'NMF')
        myObj.ROI_selection_CNN_NMF()
    end
end

%% Extract traces
% extract fluorescence traces from the raw .sbx files. Perform neuropil
% correction and calculated the dF/F (change in fluorescence with time that
% is our readout of neural activity)

myObj.sbxExtractTraces()

%% Run GLM
% run Poisson Generalized Linear Model to look at how different behavioral and task variables affects
% neural activity - using glmnet by Jerome Friedman at Stanford
% Current behavioral or task variables considered:
%     variables_to_consider = {...
%         'lick bout onset',...
%         'all other licking',...
%         'quinine',...
%         'ensure',...
%         'xyshift',...
%         'FC entire stim correct',...
%         'FC entire stim incorrect',...
%         'QC entire stim correct',...
%         'QC entire stim incorrect',...
%         'NC entire stim correct',...
%         'NC entire stim incorrect',...
%         'running',...
%         }

if strcmp(myObj.PreProcessingParameters.ROI_algorithm,'PCA/ICA')
    signals_type = 'signalsica';
elseif strcmp(myObj.PreProcessingParameters.ROI_algorithm,'NMF')
    signals_type = 'signalsnmf';
end
activity_type = 'deconvolved';
GLM_parameter_file = 'D:\Analysis_scripts\Dropbox\AndermannLab\users\ramesh\2p_preprocess\+PPPack\+hf\glm_parameters.txt';
myObj.run_glm(signals_type,activity_type,GLM_parameter_file);



