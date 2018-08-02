function write_glm_textfile_variables_used(obj,variables_to_consider,prepost_time_windows_taskvar)

% first lets define the file to write too and if it does just repoint to it
% so that don't need to recreate
runs_to_save = [];
for i = 1:length(obj.Runs_to_use)
    runs_to_save = [runs_to_save num2str(obj.Runs_to_use(i))];
end
filename_use = sprintf('%s%s_%s_runs%s_glm_variables.txt',obj.Dirs.date_mouse,obj.MouseNames,obj.Dates,runs_to_save);
if exist(filename_use)
    fileID = fopen(filename_use,'a');
else
    fileID = fopen(filename_use,'w');
end

% title
fprintf(fileID,'GLM variables used: seconds_pre, seconds_post \r\n\r\n');

% iterate through each task variable and write it
for i = 1:length(variables_to_consider)
    fprintf(fileID,'%s. %s: %s,%s \r\n',num2str(i),variables_to_consider{i},...
        num2str(prepost_time_windows_taskvar(i,1)),num2str(prepost_time_windows_taskvar(i,2)));
end

fclose(fileID);
