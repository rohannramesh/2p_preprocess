function prepost_time_windows_taskvar = get_prepost_time_windows_glm(variables_to_consider)
% get default values for time to extend pre (first value) and post (second
% value) in SECONDS
prepost_time_windows_taskvar = nan(length(variables_to_consider),2);
for i = 1:length(variables_to_consider)
    curr_task_variable = variables_to_consider{i};
    if strcmp('lick bout onset',curr_task_variable)
        prepost_time_windows_taskvar(i,:) = [1 2];
    elseif strcmp('all other licking',curr_task_variable)
        prepost_time_windows_taskvar(i,:) = [0 0];
    elseif strcmp('quinine',curr_task_variable)
        prepost_time_windows_taskvar(i,:) = [0 2];
    elseif strcmp('ensure',curr_task_variable)
        prepost_time_windows_taskvar(i,:) = [0 2];
    elseif strcmp('xyshift',curr_task_variable)
        prepost_time_windows_taskvar(i,:) = [0 0];
    elseif ~isempty(strfind(curr_task_variable,'entire'))
        prepost_time_windows_taskvar(i,:) = [0 1];
    elseif ~isempty(strfind(curr_task_variable,'onsets'))
        prepost_time_windows_taskvar(i,:) = [0 3];
    elseif ~isempty(strfind(curr_task_variable,'offsets'))
        prepost_time_windows_taskvar(i,:) = [0 3];
    elseif strcmp('running',curr_task_variable)   
        prepost_time_windows_taskvar(i,:) = [0 0];
    end
end


end