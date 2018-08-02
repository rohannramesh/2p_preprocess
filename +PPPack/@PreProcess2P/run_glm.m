function obj = run_glm(obj,signals_type,activity_type,GLM_parameters_path)
% run Poisson Generalized Linear Model to look at how different behavioral and task variables affects
% neural activity - using glmnet by Jerome Friedman at Stanford (see
% glmnet.m - Friedman et al., 2010)
% current variables in place to care about: licking, ensure, quinine,
% xyshift, all stimulus offsets, FC/QC/NC offsets, FC/QC/NC onsets, all
% onsets, all entire stim, FC/QC/NC entire stimulus, FC/QC/NC entire
% stimulus correct vs. incorrect, running, eye (see variables_to_consider_below)
%  - to adjust windows pre and post each onset go look at PPPack.hf.get_basis_set_for_GLM (bottom)
% Currently using elastic net regularization consisting of 90% L2 and 10% L1 methods
% GLM evaluated using deviance explained (Driscoll et al., 2017) and if
% want to do bootstrap test adjust tag below
% Currently training on first 75% of each run and testing on last 25%

if nargin < 4
    GLM_parameters_path = [];
end

% either pass a text file path with the number of seconds to use pre or
% post reversal or use defaults
% for formatting see glm_parameters.txt in 
if ~isempty(GLM_parameters_path)
    [variables_to_consider,prepost_time_windows_taskvar] = ...
            PPPack.hf.parse_GLM_paramters_txtfile(GLM_parameters_path);
else % defaults
    variables_to_consider = {...
        'lick bout onset',...
        'all other licking',...
        'quinine',...
        'ensure',...
        'xyshift',...
        'FC entire stim correct',...
        'FC entire stim incorrect',...
        'QC entire stim correct',...
        'QC entire stim incorrect',...
        'NC entire stim correct',...
        'NC entire stim incorrect',...
        'running',...
        }
    % this function defines the time pre and post all of the task variables to
    % extend the basis functions
    prepost_time_windows_taskvar = PPPack.hf.get_prepost_time_windows_glm(variables_to_consider);
end

% variables to input
do_bootstrap_test = 0;
ds_val = 3; % how many frames to downsample
    
% lets get the task variables for all runs
task_var = PPPack.hf.get_task_variables_for_GLM(obj.Dirs, variables_to_consider, ds_val);

% lets get the signals file to get the neural activity
nRuns = length(obj.Dirs.runs);
% to concatenate across runs
all_neural_activity = [];
% to know what frames are for training vs testing
all_train_v_test = [];
for curr_run = 1:nRuns
    % this is the directory information for a specific run
    curr_dir = obj.Dirs.runs{curr_run};
    
    % load in the appropriate signals variable
    load(strrep(curr_dir.sbx,'.sbx',['.' signals_type]),'-mat')
    % load in f2p
    load(strrep(curr_dir.sbx,'.sbx','.f2p'),'-mat')
    
    
    % number of ROIs
    nROIs = length(signals.timecourse);
    run_tc = nan(nROIs-1,floor(length(signals.timecourse(1).dff)./ds_val));
    for curr_ROI = 1:nROIs-1
        curr_vec = signals.timecourse(curr_ROI).(activity_type);
        curr_vec(curr_vec < 0) = 0;
        % now to downsample
        if ds_val ~= 1
            curr_vec = PPPack.hf.downsample_vector(curr_vec,ds_val,'mean');
        end
        run_tc(curr_ROI,:) = curr_vec;
    end
    
    % build training and testing frames - train on first 75% and test on
    % last 25% - 2 will be test and 1 will be train
    nframes_train = floor(0.75*length(curr_vec));
    curr_train_v_test_vector = 2*ones(size(curr_vec));
    curr_train_v_test_vector(1:nframes_train) = 1;
    all_train_v_test = [all_train_v_test curr_train_v_test_vector'];
    
    % concatenate across runs
    all_neural_activity = [all_neural_activity run_tc];
end



% get the basis sets
[basis_set,convolution_gauss,Beta_coeff_idx] = PPPack.hf.get_basis_set_for_GLM(obj.Dirs, variables_to_consider, task_var, ...
    ds_val, prepost_time_windows_taskvar);

% now lets run the GLM
% settings for GLM
clear opts
opts.alpha = 0.1; %used to be 0.01 L1 v L2 regularization
opts.standardize = true;
options = glmnetSet(opts);
for curr_ROI = 1:nROIs-1 % to account for neuropil
    % convolve train and testing activity with same gauss
    train_activity = all_neural_activity(curr_ROI,all_train_v_test == 1);
    train_activity = conv(train_activity,convolution_gauss,'same');
    test_activity = all_neural_activity(curr_ROI,all_train_v_test == 2);    
    test_activity = conv(test_activity,convolution_gauss,'same');
    % here is the basis set
    train_basis = sparse(basis_set(:,all_train_v_test == 1))';
    test_basis = sparse(basis_set(:,all_train_v_test == 2))';
    % remove indices with nans
    indnan = find(isnan(train_activity));
    train_activity(indnan) = [];
    if sum(train_activity == 0) == length(train_activity) | sum(isnan(train_activity)) > 0
        continue
    end
    train_basis(indnan,:) = [];
    disp(['ROI ' num2str(curr_ROI)])
    tic;CVerr = cvglmnet(train_basis,train_activity,...
        'poisson',options,'deviance',[],[],false,false,true);toc
    coeff = cvglmnetCoef(CVerr);
    Predicted_Activity_test = cvglmnetPredict(CVerr,test_basis,[],'response');
%     B = cvglmnetPredict(CVerr,train_basis,[],'response');
    % deviance quantifies how well above the mean activity of a neuron is
    % your fit doing - want more positive numbers
    [deviancevar.frac, deviancevar.D_model, deviancevar.D_null] = ...
        PPPack.hf.getDeviance(test_activity, (Predicted_Activity_test), nanmean(test_activity),'Poisson');
%     curr_ROI
%     deviancevar.frac
    GLM.models{curr_ROI} = CVerr;
    GLM.deviance{curr_ROI} = deviancevar;
    % build permutation circular shift tests to see if the fit is
    % meaningful
    if do_bootstrap_test == 1
        % bootstrap to get controlled circshift
        nIter = 50; % used to be 50
        clear deviancevar_shuff Shuffled_GLM
        for curr_boot = 1:nIter
            % have the window you circshift be 15 so that randomly shifts
            % data relative to when the actual task related variables are
            % happening
            time_range = [15];
            fr_range = round(time_range*frame_2p_metadata.meta.framerate_2p);
            circshift_val = randsample(fr_range(1):length(train_activity),1);
            CVerr_boot = cvglmnet(train_basis,circshift(train_activity',circshift_val)',...
                'poisson',options,'deviance',[],[],false,false,true);
            coeff_shuff = cvglmnetCoef(CVerr_boot);
            Predicted_Activity_test_shuff = cvglmnetPredict(CVerr_boot,test_basis,[],'response');
            [deviancevar_shuff.frac, deviancevar_shuff.D_model, deviancevar_shuff.D_null] = ...
                PPPack.hf.getDeviance(circshift(test_activity',circshift_val)', ...
                (Predicted_Activity_test_shuff), nanmean(test_activity),'Poisson');
            Shuffled_GLM.models{curr_ROI}{curr_boot} = CVerr_boot;
            Shuffled_GLM.circshift_vals{curr_ROI}{curr_boot} = circshift_val;
            Shuffled_GLM.deviance{curr_ROI}{curr_boot} = deviancevar_shuff;
        end
    else
        Shuffled_GLM = [];
    end
end
% for indexing into coefficients
GLM.opts = opts;
GLM.beta_labels = Beta_coeff_idx;

% save the GLM variable
[mouse,curr_date,runs] = PPPack.hf.get_mouse_day_run_info_from_dirs(obj.Dirs.runs{1});
save([obj.Dirs.date_mouse mouse '_' curr_date '.glm'],'GLM','Shuffled_GLM','-v7.3');

% write text file for variables used in this glm
PPPack.hf.write_glm_textfile_variables_used(obj,variables_to_consider,prepost_time_windows_taskvar);

% update log file
obj.update_log_file('GLM');
end


        


