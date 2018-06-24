function save2PPPObj(obj)
% save the object you have made

% y = inputdlg(['Enter name for saving 2P object:']);
run_name_total = [];
for i = 1:length(obj.Dirs.runs)
    run_name = obj.Dirs.runs{i}.date_mouse_run;
    ind_run = strfind(run_name,'run')
    run_name_total = [run_name_total run_name(ind_run+3:end)];
end
save([obj.Dirs.date_mouse obj.MouseNames '_' obj.Dates '_runs' run_name_total '.PreProcess'],'obj','-v7.3')