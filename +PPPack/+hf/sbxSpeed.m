function speed = sbxSpeed(mouse, date, run, server)
%SBXSPEED Return the speed of a mouse running on a 3d printed wheel from
% date from the quadrature encoder
    
    if nargin < 4, server = []; end

    speed = [];
    dirs = PPPack.hf.sbxDir(mouse, date, run, [], server);
    dirs = dirs.runs{1};
    
    if ~isempty(dirs.quad)
        quadfile = load(dirs.quad);
        running = quadfile.quad_data;
    elseif ~isempty(PPPack.hf.sbxPath(mouse, date, run, 'position', 'server', server))
        quadfile = load(PPPack.hf.sbxPath(mouse, date, run, 'position', 'server', server), '-mat');
        running = quadfile.position;
    else
        return; 
    end
    
    inf = PPPack.hf.sbxInfo(dirs.sbx);
    framerate = 30.98;
    if inf.scanmode == 1, framerate = framerate/2; end

    wheel_diameter = 14; % in cm
    wheel_tabs = 44; 
    wheel_circumference = wheel_diameter*pi;
    step_size = wheel_circumference/(wheel_tabs*2);

    instantaneous_speed = zeros(length(running), 1);
    if ~isempty(instantaneous_speed)
        instantaneous_speed(2:end) = diff(running);
        instantaneous_speed(2) = 0;
        instantaneous_speed = instantaneous_speed*step_size*framerate;
    end
    
    speed = conv(instantaneous_speed, ones(ceil(framerate/4), 1)/ceil(framerate/4), 'same')';
end

