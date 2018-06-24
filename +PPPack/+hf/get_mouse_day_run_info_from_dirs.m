function [mouse,date,run] = get_mouse_day_run_info_from_dirs(dirs)

curr_name = dirs.date_mouse_run;
ind_und = find(curr_name == '_');
date = curr_name(1:ind_und(1)-1);
mouse = curr_name(ind_und(1)+1:ind_und(2)-1);
run_string = findstr(curr_name,'run');
run = str2num(curr_name(run_string+3:end));