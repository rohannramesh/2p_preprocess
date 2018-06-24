function ch = get_ch_values_nidaq(RigName)
% These are the channels on the nidaq that can be custom to each
% independent rig
% Within the Andermann lab our 2 rigs are called Starsky and Hutch (clever
% I know) - change according to your rig name and how you hardwire the set
% up
if nargin < 1
    RigName = 'Starsky';
end


if strcmp(RigName,'Starsky')
    ch.running      = 4;
    ch.ephys        = 4; 
    ch.twoP         = 2;
    ch.visstim      = 7;
    ch.eyetrack     = 7;
    ch.bodytrack    = 7; 
    ch.licking      = 6;
    ch.ensure       = 5;
    ch.OTB_visstim  = 7;
    ch.quinine      = 3;
elseif strcmp(RigName,'Hutch')
    ch.running      = 5;
%     ch.ephys        = 3; 
    ch.twoP         = 6;
    ch.visstim      = 7;
    ch.eyetrack     = 6;
    ch.bodytrack    = 8; 
%     ch.epi          = 7;
    ch.licking      = 4;
    ch.ensure       = 3;
    ch.OTB_visstim  = 2;
    ch.quinine      = 9;
end