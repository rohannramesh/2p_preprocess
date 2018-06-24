function processForJavascript(mouse, date, runs, force, axon, server, nmf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 6, server = []; end
    if nargin < 3 || isempty(runs), error('Need runs'); end    
    if nargin < 4 || isempty(force), force = false; end
    if nargin < 5 || isempty(axon), axon = false; end
    if nargin < 7 || isempty(nmf), nmf = false; end
    
    icarun = runs(end);

    icapath = PPPack.hf.sbxPath(mouse, date, icarun, 'ica', 'server', server);
    nmftext = '';
    if nmf || isempty(icapath)
        icapath = PPPack.hf.sbxPath(mouse, date, icarun, 'icanmf', 'server', server);
        nmftext = '_nmf';
    end
    
    % this is the directory for the Online ROI identification GUI
    datapath = '\\tolman\webdata\mousedata\';
    
    mdp = sprintf('%s%s_%s_%03i%s\\', datapath, mouse, date, icarun, nmftext);
    zeropos = 5;  % height above zero relative to below
    minim = 20;  % Minimum image width
    erosions = [1.0, 0.85, 0.7, 0.5, 0.35, 0.2, 0.12, 0.05];
    
    if ~exist(mdp) || force
        mkdir(mdp)
        ica = load(icapath, '-mat');
        
        mnimpath = sprintf('%smean-image.jpg', mdp);
        PPPack.hf.saveScaledImage(mnimpath, ica.premasks.movm);
        
        roipath = sprintf('%srois.txt', mdp);
        fp = fopen(roipath, 'w');
        fprintf(fp, '%i', length(ica.premasks.ica));
        fclose(fp);
        
        impos = zeros(length(ica.premasks.ica), 3);
        masks = cell(1, length(ica.premasks.ica));
        erosionmasks = cell(1, length(ica.premasks.ica));
        
        for tr = 1:length(ica.premasks.ica)
            trpath = sprintf('%strace-%04i.jpg', mdp, tr);
            impath = sprintf('%sfilter-%04i.png', mdp, tr);
            erpath = sprintf('%smask-%04i.json', mdp, tr);
            
            PPPack.hf.saveTraceImage(trpath, ica.premasks.ica(tr).trace, zeropos);
            [impos(tr, 1), impos(tr, 2), impos(tr, 3)] = ...
                PPPack.hf.saveFilterImage(impath, ica.premasks.ica(tr).filter, minim);
            [erosionmasks{tr}, finerosions] = ...
                PPPack.hf.saveErosionMasks([], ica.premasks.ica(tr).filter, erosions, ~axon);
            masks{tr} = PPPack.hf.getErosionMask(ica.premasks.ica(tr).filter, 0.7, ~axon);
        end
        
        impospath = sprintf('%simage-positions.json', mdp);
        fp = fopen(impospath, 'w');
        fprintf(fp, '{"positions":[');
        for tr = 1:length(ica.premasks.ica)
            if tr > 1, fprintf(fp, ','); end
            fprintf(fp, '[%i,%i,%i]', impos(tr, 1), impos(tr, 2), impos(tr, 3));
        end
        fprintf(fp, ']}');
        fclose(fp);
        
        ermaskpath = sprintf('%serosion-masks.json', mdp);
        fp = fopen(ermaskpath, 'w');
        fprintf(fp, '{"erosions":[');
        for (e = 1:length(finerosions))
            if (e > 1), fprintf(fp, ','); end
            fprintf(fp, '%0.2f', finerosions(e));
        end
        fprintf(fp, '], "masks":[');
        for tr = 1:length(ica.premasks.ica)  % For cell
            if tr > 1, fprintf(fp, ','); end
            fprintf(fp, '[');
            for i = 1:length(erosionmasks{tr})  % For erosion level
                if i > 1, fprintf(fp, ','); end
                fprintf(fp, '[');
                for j = 1:length(erosionmasks{tr}{i})  % Mask array
                    if j > 1, fprintf(fp, ','); end
                    fprintf(fp, '%i', erosionmasks{tr}{i}(j));
                end
                fprintf(fp, ']');
            end
            fprintf(fp, ']');
        end
        fprintf(fp, ']}');
        fclose(fp);
        
        overlappath = sprintf('%soverlaps.json', mdp);
        PPPack.hf.saveOverlapMasks(overlappath, masks);
    end
end

