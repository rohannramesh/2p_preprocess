function base = sbxScanbase(server)
%SBXSCANBASE Hard codes directories

    % Set base path depending on server
    if nargin < 1 || isempty(server) || strcmpi(PPPack.hf.hostname, server)
        if strcmp(PPPack.hf.hostname, 'Megatron')
            base = 'D:\twophoton_data\2photon\scan\';
        elseif strcmp(PPPack.hf.hostname, 'Atlas')
            base = 'E:\twophoton_data\2photon\raw\';
        elseif strcmp(PPPack.hf.hostname, 'BeastMode')
            base = 'S:\twophoton_data\2photon\scan\';
        elseif strcmp(PPPack.hf.hostname, 'Sweetness')
            base = 'D:\2p_data\scan\';
        elseif strcmpi(PPPack.hf.hostname, 'santiago')
            base = 'D:\2p_data\scan\';
        end
    else
        if strcmpi(server, 'santiago')
            base = '\\santiago\2p_data\scan\';
        elseif strcmpi(server, 'sweetness')
            base = '\\sweetness\2p_data\scan\';
        elseif strcmpi(server, 'megatron')
            base = '\\megatron\2photon\scan\';
        elseif strcmpi(server, 'storage') && strcmp(hostname, 'Megatron')
            base = 'E:\scan\';
        elseif strcmpi(server, 'storage')
            base = '\\megatron\E\scan\';
        elseif strcmpi(server, 'Anastasia')
            base = get_anastasia_pathname(PPPack.hf.username);
        end
    end
end

function base = get_anastasia_pathname(curr_username)
    % to get the base directory for anastasia which is organized
    % differently
    anastasia_base = '\\ANASTASIA\data\2p\';
    if strcmpi(curr_username,'kmcguir2') || strcmpi(curr_username,'kmcguire')
        curr_name = 'kelly\';
    elseif strcmpi(curr_username,'cburgess')
        curr_name = 'kelly\';
    elseif strcmpi(curr_username,'rramesh1') || strcmpi(curr_username,'rramesh')
        curr_name = 'kelly\';
    elseif strcmpi(curr_username,'cburgess')
        curr_name = 'kelly\';
    end
    base = [anastasia_base curr_name];
end


