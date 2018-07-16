classdef PreProcess2P < handle
% This class will perform all preprocessing steps for twophoton microscopy
% data collected using the Scanbox system and file structure (.sbx files).
% All data is expected to be in the format:
% Mouse folder name (e.g. RR6), Date of collection (180101_RR6), Individual
% runs (180101_RR6_run1) with the last folder containing all of the data
% created withing the scanbox system (.sbx and info file) provided by
% Neurolabware. If you do not have this file organization edit here: PPPack.hf.sbxScanbase
% Pre-processing of the data includes movie registration, tiff creation to
% check for proper registration, ROI identification, GUIs / CNN for
% selection of ROIs, trace extraction, and GLM
% (assumes target run is the first one for registrations)
%
%
% INPUTS (two options):
% 1. Empty square brackets []: a gui will pop up and ask you to enter the
% MouseName, the Date, and the runs to analyze
% 2. Directly pass MouseName, the Date, and the runs to analyze



    properties
        MouseNames                  % the name of the mouse to be analyzed
        Dates                       % the date to analyze
        Runs_to_use                 % which runs to analyze
        Target_run                  % the run to register all runs on that day to
        Dirs                        % the directories for all of the data with defined pathnames for all variables of interest
        PreProcessingParameters     % preprocessing parameters include rig imaged, registration type, ROI type (neuron vs. axon), ROI selection algorithm, Neuropil correction type 
    end
    
    methods
        sbxRegister(obj);
        sbxTestRegistrationImageStack(obj)
        sbxframe2pmetadata(obj)
        sbxTrialInfo(obj)
        PCAICA_ROI_identification(obj)
        CNMF_ROI_identification(obj)
        run_ROI_selection_gui(obj)
        sbxExtractTraces(obj)
        ROI_selection_CNN_NMF(obj)
        ROI_selection_CNN_PCAICA(obj)
        run_glm(obj,signals_type,activity_type)
        save2PPPObj(obj)
        update_log_file(obj,tag)
        function obj = PreProcess2P(varargin)
            %% Empty input --> will populate via an inputdlg
            if isempty(varargin{1}) || isstruct(varargin{1})
                prompt = {'Enter mouse name:','Enter date to analyze:', ...
                    'Enter runs to analyze:'};
                title = 'Input';
                dims = [1 35];
                definput = {'RR1','180101','1:4'};
                answer = inputdlg(prompt,title,dims,definput);
                obj.MouseNames = answer{1};
                obj.Dates = answer{2};
                obj.Runs_to_use = str2num(answer{3});
                if length(obj.Runs_to_use) > 6
                    warning('Might struggle with this many runs at one time')
                end
                obj.Target_run = obj.Runs_to_use(1);
                % now have a function that compares any specific preprocessing
                % formatting requests with a default expectation based off
                % username and standard lab default practices
                obj.PreProcessingParameters = PPPack.hf.get_PreProcessingParameters(varargin{1});
            else % if pass MouseName, Date, and Runs
                obj.MouseNames = varargin{1};
                obj.Dates = varargin{2};
                obj.Runs_to_use = varargin{3};
                obj.Target_run = obj.Runs_to_use(1);
                % now have a function that compares any specific preprocessing
                % formatting requests with a default expectation based off
                % username and standard lab default practices
                if length(varargin) == 4
                    obj.PreProcessingParameters = PPPack.hf.get_PreProcessingParameters(varargin{4});
                else
                    obj.PreProcessingParameters = PPPack.hf.get_PreProcessingParameters([]);
                end
            end
            % now have to build the directories
            obj.Dirs = PPPack.hf.sbxDir(obj.MouseNames,obj.Dates,obj.Runs_to_use,...
                obj.Target_run,obj.PreProcessingParameters.server);
        end
         
    end
end