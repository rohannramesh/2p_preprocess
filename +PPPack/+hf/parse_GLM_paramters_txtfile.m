function [variables_to_consider,prepost_time_windows_taskvar] = ...
    parse_GLM_paramters_txtfile(parameter_path)

fileID = fopen(parameter_path,'r');
   
InputText = textscan(fileID,'%s',1,'delimiter','\n');  % Read header lines
Labels_prepost = textscan(fileID,'%s%s%s',1,'delimiter',':,'); % labels
% disp(InputText{1});                        % Display header lines

% read each data row
dataRow = 1;
idx = 1;
clear variables_to_consider prepost_time_windows_taskvar
while dataRow
  InputText = textscan(fileID,'%s%f%f', 1,...    % Read data row
      'delimiter',':,');  
  if findstr(InputText{1}{1},'*') % marking the end of block to use with *EOB
      dataRow = 0;
  else
      variables_to_consider{idx} = InputText{1}{1};
      prepost_time_windows_taskvar(idx,1) = InputText{2};
      prepost_time_windows_taskvar(idx,2) = InputText{3};
      idx = idx + 1;
  end
end

fclose(fileID);