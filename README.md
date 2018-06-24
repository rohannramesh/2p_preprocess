# 2p_preprocess

This is a Pre-Processing package for 2-photon imaging data sets. Here, I will describe briefly the approach taken at each step, the parameters that should be set depending on your particular imaging setup, and decisions that need to be made in terms of how you want to preprocess your data. Within the PPPack (Pre-Process Package) there are three main folders: the PreProcess2P class, the helper functions package (hf), the constrained non-negative matrix factorization (NMF) package, and the PCA/ICA package. These last two folders contain the functions and pre-trained convolutional neural networks needed to identify putative regions of interest (ROIs) and for ROI selection


## Setup, Prerequisites, and Installing

These instructions will get you a copy of the project up and running on your local machine in Matlab.

All code has been tested and optimized in Matlab 2015b+. It is unclear if
all functions will work the same in earlier versions. Code was originally implemented on Windows Servers.

Dependencies: Statistics and Machine Learning Toolbox, Image Processing toolbox.
Recommended but not necessary: Parallel compution toolbox (for large data sets) and Signal processing toolbox 

Install python v 3.6 and NumPy, Keras, and Tensorflow if you would like to use the pre-trained convolutional neural networks for ROI selection and make sure Matlab acknowledges Python's existence using pyversion. 

From within +PPPack move ca_source_extraction_master and glmnet_matlab into your path. Also add both ROI selection GUI's for both PCA/ICA and NMF to your path (both .fig and m file). All folders with a + in front of it are protected and can not be added to the path. 
For deconvolution, we are using an algorithm that requires the CVX library, download here:  http://cvxr.com/cvx/download/

### Preprocessing pipeline

PreProcessing Steps - to run look at Run_PreProcess_Pipeline.m.
Additional details also included about each step included in this script
1. Define your PreProcessing object
2. Save object
3. Register data
4. Build frame metadata and trial information - (optional) if using MonkeyLogic behavioral training setup (Asaad and Esakandar, 2008)
5. ROI identification - PCA/ICA (Mukamel and Schniter, 2009) or CNMF (Pnevmatikakis et al., 2016)
6. ROI extraction - pre-trained convolutional neural network or ROI selection GUI (for semi-automated manual selection)
7. Signal fluorescence trace extraction (neuropil correction and dFF calculation)
8. Run Generalized Linear Model (optional)


## Imaging Protocol and File Structure

Each 2P imaging session can contain multiple imaging runs. In the Andermann lab we collect ~4 30 min imaging runs at a single plane. A separate object should be created for each new plane imaged. Within a given imaging plane there will be multiple active neurons. The goal of this preprocessing pipeline is to register all imaging runs to a common target with subpixel resolution, to identify those neurons that are active during an imaging run, and to extract the activity of these neurons for all imaging runs on a given day. Our lab uses the relatively common Scanbox imaging setup and so within the folder for an imaging run there will be the following files:
1. .sbx - raw data from that imaging run (still have to make flexible
          for tiff setup)
2. .ephys - nidaq data if aligning animal behavior to 2P frames
3. .mat - info file for 2P imaging run (paired with .sbx file)
4. eye.mat - file with eye tracking data
5. .quadrature - information about animal's running

The pre-processing files will build the following files (load each filetype as a common .mat file- format used is for easy file identification)
1. .align - file containing the subpixel shifts for each plane to 
            register each run to a common target
2. .f2p and .TrialVar - file for information about stimulus and behavior on each frame
3. .ica or .nmf - output from PCA/ICA or CNMF ROI identification algorithm respectively
4. .icamasks or .nmfmasks - output after ROI selection process (GUI or CNN)
5. .npilmasks - masks for each neuropil mask for each ROI identified
6. .signalsica or .signalsnmf - signals variable containing the extracted timecourses

Folder organization we used:
All of the data for each run was saved in a single folder, with multiple folders for multiple runs for each day.
An example path to a single run folder:
'S:\twophoton_data\2photon\scan\OA27\170524_OA27\170524_OA27_run4.' 
The important considerations are Base\MouseName\Date\Runs_to_analyze. We define base in PPPack.hf.sbxScanBase and everything else in PPPack.hf.sbxDir From personal experience we recommend this structure, but feel free to adjust accordingly.


## Authors

* Rohan Ramesh - PhD. student, Andermann lab
