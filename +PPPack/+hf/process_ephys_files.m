function nidaq = process_ephys_files(pathname, varargin)

if length(varargin) == 2
    TwoP_framerate = varargin{1};
    nchannels = varargin{2};
else
    % load in global variable info
    sbxpath = [pathname(1:end-6) '.sbx'];
    info = PPPack.hf.sbxInfo(sbxpath);
    if info.scanmode == 1
        TwoP_framerate = 15.5;
    elseif info.scanmode == 0
        TwoP_framerate = 31;
    end
end
% first lets load in variable
    A = fopen(pathname);
    B = fread(A,'float');
% get the number of channels
    nidaq.Fs = 2000;
    display(['Assuming Fs is 2000']); % this is based off the Scanbox system and its common implementation
if nargin == 1
    convert_to_sec = length(B)./nidaq.Fs;
    length_mov = (info.max_idx+1)./TwoP_framerate;
    nchannels = round(convert_to_sec./length_mov);
end

% this will put the ehpys or nidaq file into an output where each column of
% the variable data corresponds to one channel recording pulses
display(['Assuming you have ' num2str(nchannels) ' channels and unshuffling accordingly']);
nidaq.data = PPPack.hf.unshuffle_array(B',nchannels);

% make the timeStamps variable artificially
nidaq.timeStamps = ((1:size(nidaq.data,1))-1)'./nidaq.Fs;

% Close the file.
fclose(A);
