function obj = update_log_file(obj,tag)
% Create a running description of what has been run for each step of
% processing so that you can easily open a text file and see what has been
% done so far. This will save as MouseName_Date_runsXYZ_log.txt in the
% folder of a given day (Date_MouseName)

% first lets define the file to write too and if it does just repoint to it
% so that don't need to recreate
runs_to_save = [];
for i = 1:length(obj.Runs_to_use)
    runs_to_save = [runs_to_save num2str(obj.Runs_to_use(i))];
end
filename_use = sprintf('%s%s_%s_runs%s_log.txt',obj.Dirs.date_mouse,obj.MouseNames,obj.Dates,runs_to_save);
if exist(filename_use)
    fileID = fopen(filename_use,'a');
else
    fileID = fopen(filename_use,'w');
end

% these are the tags that are called in their respective functions so that
% a user can see what has been run for a given run on a given day
% if build new functions in the class then add here
if strcmp(tag,'registration')
    fprintf(fileID,'Registration: %s runs %s\r\n',obj.PreProcessingParameters.registration,runs_to_save);
elseif strcmp(tag,'frame_2p_metadata')
    fprintf(fileID,'f2p created runs %s\r\n',runs_to_save);
elseif strcmp(tag,'trial_info')
    fprintf(fileID,'TrialVar created runs %s\r\n',runs_to_save);
elseif strcmp(tag,'PCAICA_ROI_identification')
    fprintf(fileID,'.ica created runs %s\r\n',runs_to_save);
elseif strcmp(tag,'NMF_ROI_identification')
    fprintf(fileID,'.nmf created runs %s\r\n',runs_to_save);
elseif strcmp(tag,'CNN_PCAICA')
    fprintf(fileID,'CNN run for .ica runs %s\r\n',runs_to_save);    
elseif strcmp(tag,'CNN_NMF')
    fprintf(fileID,'CNN run for .nmf runs %s\r\n',runs_to_save);  
elseif strcmp(tag,'ROI_selection_GUI')
    fprintf(fileID,'GUI run for %s runs %s\r\n',obj.PreProcessingParameters.ROI_algorithm,runs_to_save);
elseif strcmp(tag,'trace_extraction')
    fprintf(fileID,'traces extracted for runs %s\r\n',runs_to_save);  
elseif strcmp(tag,'GLM')
    fprintf(fileID,'GLM run for runs %s\r\n',runs_to_save);  
else
    sprintf('Did not recognize tag')
end



fclose(fileID);
