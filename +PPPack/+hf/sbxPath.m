function path = sbxPath(mouse, date, run, type, varargin)
%SBXPATH Gets the path to most types of sbx files. Requires run.

    p = inputParser;
    addOptional(p, 'estimate', false);  % Whether to give an estimate of where the path would be if it does not exist
    addOptional(p, 'server', []);  % Which server to call, if not the current server
    addOptional(p, 'pmt', []);  % Which PMT, in the case of sbxreg and sbxclean
    parse(p, varargin{:});
    p = p.Results;
    
    path = [];
    
    % If date is an integer, convert to string
    if ~ischar(date), date = num2str(date); end

    % Load files from scanbox directories
    dirsf = PPPack.hf.sbxDir(mouse, date, run, [], p.server);
    dirs = dirsf.runs{1};
    
    searchdir = dirs.path;
    if isempty(searchdir), return; end
    if isempty(dir(searchdir)), return; end
    
    % Get path
    switch type
        case 'sbx'
            path = dirs.sbx;
        case 'clean'
            if ~isempty(dirs.sbx)
                [a, b, ~] = fileparts(dirs.sbx);
                cpath = [a '\' b '.sbxclean'];
                if exist(cpath, 'file'), path = cpath; end
            end
        case 'info'
            if ~isempty(dirs.sbx), path = dirs.sbx(1:end - 4); end
        case 'stim'
            path = [dirs.path dirs.sbx_name '.f2p'];
        case 'simpcell'
            fname = sprintf('%s_%s_%03i.simpcell', mouse, date, run);
            path = [dirsf.date_mouse fname];
        case 'simpglm'
            fname = sprintf('%s_%s.simpglm', mouse, date);
            path = [dirsf.date_mouse fname];
        case 'bhv'
            path = dirs.ml;
        case 'ephys'
            path = dirs.nidaq;
        case 'pupil'
            if ~isempty(dirs.sbx)
                checkpath = [dirs.sbx(1:end - 4) '_eye.mat']; 
                if exist(checkpath, 'file'), path = checkpath; end
            end
        case 'quad'
            path = dirs.quad;
            
        case 'xyreg'
            type = 'sbxreg';
            fs = dir(searchdir);
            for i=1:length(fs)
                [~, fname, ext] = fileparts(fs(i).name);
                if strcmp(ext, sprintf('.%s', type))
                    if isempty(p.pmt) || ~isempty(strfind([fname ext], sprintf('_xyreg-%i.sbxreg', p.pmt)))
                       path = sprintf('%s%s', searchdir, fs(i).name);
                    end
                end
            end
        case 'demonsreg'
            type = 'sbxreg';
            fs = dir(searchdir);
            for i=1:length(fs)
                [~, fname, ext] = fileparts(fs(i).name);
                if strcmp(ext, sprintf('.%s', type))
                    if isempty(p.pmt) || ~isempty(strfind([fname ext], sprintf('_demonsreg-%i.sbxreg', p.pmt)))
                        path = sprintf('%s%s', searchdir, fs(i).name);
                    end
                end
            end
        case 'xyregclean'
            type = 'sbxclean';
            fs = dir(searchdir);
            for i=1:length(fs)
                [~, fname, ext] = fileparts(fs(i).name);
                if strcmp(ext, sprintf('.%s', type))
                    if isempty(p.pmt) || ~isempty(strfind(fname, sprintf('_xyreg-%i', p.pmt)))
                        path = sprintf('%s%s', searchdir, fs(i).name);
                    end
                end
            end
        otherwise
            fs = dir(searchdir);
            [fs,~]=nestedSortStruct(fs,'date'); %AL added 180322
            
            for i=1:length(fs)
                [~, fname, ext] = fileparts(fs(i).name);
                if strcmp(ext, sprintf('.%s', type))
                    %if isempty(p.pmt) || ~isempty(strfind(fname, sprintf('_reg-%i', p.pmt)))
                    if isempty(p.pmt) || ~isempty(strfind([fname '.' type], sprintf(['_reg-%i.' type], p.pmt))) % AL changed 180322
                       path = sprintf('%s%s', searchdir, fs(i).name);
                    end
                end
            end
    end
    
    if nargin > 4 && ((islogical(p.estimate) && p.estimate) || (~islogical(p.estimate) && ~isempty(p.estimate)))
        path = [dirs.base '.' type];
    end

end

