function [out, varargout] = sbxLoad(mouse, date, run, type, server, varargin)
% SBXLOAD loads any type of important file of the scanbox format. Has
% minimal path assumptions, but depends on sbxDir.
% if enter in sbx then loading the sbx file and then specify 'subpixel' or
% 'whole pixel' as first entry of varargin

    out = [];

    if nargin < 4
        disp('ERROR: Call with mouse, date, run, and type of file.');
        return
    elseif nargin < 5
        server = [];
    end
    
    % If date is an integer, convert to string
    if ~ischar(date), date = num2str(date); end
    
    % Get the file path
    path = PPPack.hf.sbxPath(mouse, date, run, type, 'server', server);
    
    % Share the results
    if isempty(path)
        disp('WARNING: File not found');
        return;
    else
        %disp(sprintf('File type %s found at %s.', type, path));
    end
    
    dirsf = sbxDir(mouse, date, run, [], server);
    dirs = dirsf.runs{1};
    
    % Read the necessary filetypes differently
    switch type
        case 'sbx'
            out = PPPack.hf.sbxLoadRegRun(dirs,varargin{1});
        case 'info'
            out = PPPack.hf.sbxInfo(path);
        case 'stim_master'
            out = load(strrep(dirs.sbx,'.sbx','.TrialVar'));
        case 'bhv'
            out = PPPack.hf.bhv_read(dirs.ml);
        case 'ephys'
            out = PPPack.hf.process_ephys_files(strrep(dirs.sbx,'.sbx','.ephys'));
        otherwise
            if exist(path, 'file')
                out = load(path, '-mat');
            end
    end
end