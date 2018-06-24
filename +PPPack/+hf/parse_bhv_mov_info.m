function stim = parse_bhv_mov_info(stim,curr_Timing_Script,count)
% This function will be used to parse the bhv file and extract information
% about stimulus orientation etc. shown in each movied presented to the
% animal


% re-derive Condition and unique condition
Condition = stim.ConditionNumber; % Vector that lists which condition occurred each trial
unique_Condition = unique(Condition);
unique_TaskObject = stim.TaskObject(unique_Condition,:);
% these are the images/ movies for this condition
% if blanks this won't work bc no blank movie so have to hardcode that part
if strfind(curr_Timing_Script{1},'Blank')
    % need to format this way bc of later steps - 360 degree stimulus
    % will be used to denote a blank stimulus
    curr_task_objects = '_360_';
else
    for curr_task_object_ind = 1:size(unique_TaskObject,2)
        curr_object = unique_TaskObject{count,curr_task_object_ind};
        % lets identify if movie or a sound file
        if isempty(curr_object)
            continue
        elseif strfind(curr_object,'snd')
            curr_task_objects = '_-1_';
        elseif isempty(strfind(curr_object,'Mov'))
            continue
        else
            % only look in parentheses
            ind_par = find(curr_object == '(');
            curr_task_objects = curr_object(ind_par+1:end-1);
            % throw out everything after last comma bc this is specifying
            % location on screen
            curr_task_objects(find(curr_task_objects == ',',1,'first')+1:end) = [];
        end
    end
end

ind_CS = find(Condition == unique_Condition(count));
% now just put those things that don't change - these are constants across
% all movies
stim.frame.sf(ind_CS)          = 0.04;
stim.frame.tf(ind_CS)          = 2;
stim.frame.location(ind_CS)    = -1;
if strfind(curr_task_objects,'Contr') % if have specified a contrast in the movie will draw the value out
    ind_contr = strfind(curr_task_objects,'Contr_');
    ind_comma = find(curr_task_objects == ',');
    stim.frame.contrast(ind_CS)    = str2num(curr_task_objects(ind_contr+6:ind_comma(1)-1));
else
    stim.frame.contrast(ind_CS)    = -0.8; % otherwise assuming contrast of 0.8 and negative implies square vs sine gratings
end
% now for CS ori
ori_possibilities = [-1 0 45 90 135 180 225 270 315 360]; % -1 for auditory and others are orientation possibilities
clear foundStr
try for i = 1:length(ori_possibilities)  
        % if one of the possible oris are part of the string identify it bc
        % offset by underscores or last one is a comma
        % most common format is Mov_0_Contr_0.1.avi
        tmpRU = strfind(curr_task_objects,['_' num2str(ori_possibilities(i)) '_']);
        tmpRC = strfind(curr_task_objects,['_' num2str(ori_possibilities(i)) ',']);
        tmpR = [tmpRU tmpRC];
        if tmpR
            stim.frame.orientation(ind_CS) = ori_possibilities(i);
            foundStr = 1;
        end
        % if you have gone through all possibilities and there is no
        % matching string then artificially make an error and hop to the
        % hardcoded options
        if i == length(ori_possibilities) && ~exist('foundStr')
            falseError
        end
    end
catch err
    if ~isempty(findstr(curr_Timing_Script{1},'CSp'))
        stim.frame.orientation(ind_CS) = 0;
            warning('Assuming hardcoded ori')
    elseif ~isempty(findstr(curr_Timing_Script{1},'CSm'));
        stim.frame.orientation(ind_CS) = 270;
            warning('Assuming hardcoded ori')
    elseif ~isempty(findstr(curr_Timing_Script{1},'CSn'));
        stim.frame.orientation(ind_CS) = 135;
            warning('Assuming hardcoded ori')
    elseif ~isempty(findstr(curr_Timing_Script{1},'Blank'));
        stim.frame.orientation(ind_CS) = 360;
            warning('Assuming hardcoded ori')
    end
end



