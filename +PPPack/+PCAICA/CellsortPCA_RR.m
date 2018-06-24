function [mixedsig, mixedfilters, CovEvals, covtrace, movm, ...
    movtm] = CellsortPCA_RR(mov, flims, nPCs, dsamp, outputdir, badframes)
% [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(fn, flims, nPCs, dsamp, outputdir, badframes)
%
% CELLSORT
% Read TIFF movie data and perform singular-value decomposition (SVD)
% dimensional reduction.
%
% Inputs:
%   fn - movie file name. Must be in TIFF format.
%   flims - 2-element vector specifying the endpoints of the range of
%   frames to be analyzed. If empty, default is to analyze all movie
%   frames.
%   nPCs - number of principal components to be returned
%   dsamp - optional downsampling factor. If scalar, specifies temporal
%   downsampling factor. If two-element vector, entries specify temporal
%   and spatial downsampling, respectively.
%   outputdir - directory in which to store output .mat files
%   badframes - optional list of indices of movie frames to be excluded
%   from analysis
%
% Outputs:
%   mixedsig - N x T matrix of N temporal signal mixtures sampled at T
%   points.
%   mixedfilters - N x X x Y array of N spatial signal mixtures sampled at
%   X x Y spatial points.
%   CovEvals - largest eigenvalues of the covariance matrix
%   covtrace - trace of covariance matrix, corresponding to the sum of all
%   eigenvalues (not just the largest few)
%   movm - average of all movie time frames at each pixel
%   movtm - average of all movie pixels at each time frame, after
%   normalizing each pixel deltaF/F
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

tic
% fprintf('-------------- CellsortPCA %s: %s -------------- \n', date, fn)

%-----------------------


%-----------------------
% Check inputs
% if ~exist(
%     error('Invalid input file name.')
% end
if (nargin<2)||(isempty(flims))
%     nt_full = tiff_frames(fn);
    nt_full = size(mov,3);
    flims = [1,nt_full];
end

useframes = setdiff((flims(1):flims(2)), badframes);
nt = length(useframes);

if nargin<3 || isempty(nPCs)
    nPCs = min(150, nt);
end
if nargin<4 || isempty(dsamp)
    dsamp = [1,1];
end
if nargin<5 || isempty(outputdir)
    outputdir = [pwd,'/cellsort_preprocessed_data/'];
end
if nargin<6
    badframes = [];
end
if isempty(dir(outputdir))
    mkdir(pwd, '/cellsort_preprocessed_data/')
end
if outputdir(end)~='/';
    outputdir = [outputdir, '/'];
end

% Downsampling
if length(dsamp)==1
    dsamp_time = dsamp(1);
    dsamp_space = 1;
else
    dsamp_time = dsamp(1);
    dsamp_space = dsamp(2); % Spatial downsample
end


% Get movie info
[pixw,pixh,~] = size(mov);

npix = pixw*pixh;

fprintf('   %d pixels x %d time frames;', npix, nt)
if nt<npix
    fprintf(' using temporal covariance matrix.\n')
else
    fprintf(' using spatial covariance matrix.\n')
end

% Create covariance matrix
if nt < npix
    [covmat, mov, movm, movtm] = create_tcov(mov, pixw, pixh, useframes, nt, dsamp);
else
    [covmat, mov, movm, movtm] = create_xcov(mov, pixw, pixh, useframes, nt, dsamp);
end
covtrace = trace(covmat) / npix;
movm = reshape(movm, pixw, pixh);

if nt < npix
    % Perform SVD on temporal covariance
    [mixedsig, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix);

    % Load the other set of principal components
    [mixedfilters] = reload_moviedata(pixw*pixh, mov, mixedsig, CovEvals);
else
    % Perform SVD on spatial components
    [mixedfilters, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix);

    % Load the other set of principal components
    [mixedsig] = reload_moviedata(nt, mov', mixedfilters, CovEvals);
end
mixedfilters = reshape(mixedfilters, pixw,pixh,nPCs);

% firstframe_full = imread(fn,'tiff',1,'Info',info); %MJLM
firstframe_full = mov(:,:,1); %MJLM


firstframe = firstframe_full;
if dsamp_space>1
    firstframe = imresize(firstframe, size(mov(:,:,1)),'bilinear');
end

%------------
% Save the output data
% save(fnmat,'mixedfilters','CovEvals','mixedsig', ...
%     'movm','movtm','covtrace')
fprintf(' CellsortPCA: saving data and exiting; ')
toc

    function [covmat, mov, movm, movtm] = create_xcov(mov, pixw, pixh, useframes, nt, dsamp)
        %-----------------------
        % Load movie data to compute the spatial covariance matrix

        npix = pixw*pixh;
        
        % Downsampling
        if length(dsamp)==1
            dsamp_time = dsamp(1);
            dsamp_space = 1;
        else
            dsamp_time = dsamp(1);
            dsamp_space = dsamp(2); % Spatial downsample
        end

%         if (dsamp_space==1)
%             mov = zeros(pixw, pixh, nt);
% 			
%             for jjind=1:length(useframes)
%                 jj = useframes(jjind);
%                 %MJLM speedup:
% 				%mov(:,:,jjind) = imread(fn,jj);
% (				mov(:,:,jjind) = imread(fn,'tiff',jj,'Info',info);
% 				
%                 if mod(jjind,500)==1
%                     fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
%                     toc
%                 end
%             end
%         else
%             [pixw_dsamp,pixh_dsamp] = size(imresize( imread(fn,'tiff',1,'Info',info), 1/dsamp_space, 'bilinear' ));
%             mov = zeros(pixw_dsamp, pixh_dsamp, nt);
%             for jjind=1:length(useframes)
%                 jj = useframes(jjind);
%                 mov(:,:,jjind) = imresize( imread(fn,'tiff',jj,'Info',info), 1/dsamp_space, 'bilinear' );
%                 if mod(jjind,500)==1
%                     fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
%                     toc
%                 end
%             end
%         end

        fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
        toc
        mov = reshape(mov, npix, nt);

        % DFoF normalization of each pixel
        movm = mean(mov,2); % Average over time
        movmzero = (movm==0);
        movm(movmzero) = 1;
        mov = mov ./ (movm * ones(1,nt)) - 1; % Compute Delta F/F
        mov(movmzero, :) = 0;

        if dsamp_time>1
            mov = filter(ones(dsamp,1)/dsamp, 1, mov, [], 2);
            mov = downsample(mov', dsamp)';
        end

        movtm = mean(mov,2); % Average over space
        clear movmzeros

        c1 = (mov*mov')/size(mov,2);
        toc
        covmat = c1 - movtm*movtm';
        clear c1
    end

    function [covmat, mov, movm, movtm] = create_tcov(mov, pixw, pixh, useframes, nt, dsamp)
        %-----------------------
        % Load movie data to compute the temporal covariance matrix
        npix = pixw*pixh;

        % Downsampling
        if length(dsamp)==1
            dsamp_time = dsamp(1);
            dsamp_space = 1;
        else
            dsamp_time = dsamp(1);
            dsamp_space = dsamp(2); % Spatial downsample
        end


%         fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
        toc
        mov = reshape(mov, npix, nt);

        % DFoF normalization of each pixel
        movm = mean(mov,2); % Average over time
        movmzero = (movm==0); % Avoid dividing by zero
        movm(movmzero) = 1;
        mov = mov ./ (movm * ones(1,nt)) - 1;
        mov(movmzero, :) = 0;

        if dsamp_time>1
            mov = filter(ones(dsamp,1)/dsamp, 1, mov, [], 2);
			% MJLM: Changed this to downsampling with averaging:
            %mov = downsample(mov', dsamp)';
			mov = downsampleWithAvg(mov, dsamp);
        end

        c1 = (mov'*mov)/npix;
        movtm = mean(mov,1); % Average over space
        covmat = c1 - movtm'*movtm;
        clear c1 
    end

    function [mixedsig, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix)
        %-----------------------
        % Perform SVD

        covtrace = trace(covmat) / npix;

        opts.disp = 0;
        opts.issym = 'true';
        if nPCs<size(covmat,1)
            [mixedsig, CovEvals] = eigs(covmat, nPCs, 'LM', opts);  % pca_mixedsig are the temporal signals, mixedsig
        else
            [mixedsig, CovEvals] = eig(covmat);
            CovEvals = diag( sort(diag(CovEvals), 'descend'));
            nPCs = size(CovEvals,1);
        end
        CovEvals = diag(CovEvals);
        if nnz(CovEvals<=0)
            nPCs = nPCs - nnz(CovEvals<=0);
            fprintf(['Throwing out ',num2str(nnz(CovEvals<0)),' negative eigenvalues; new # of PCs = ',num2str(nPCs),'. \n']);
            mixedsig = mixedsig(:,CovEvals>0);
            CovEvals = CovEvals(CovEvals>0);
        end

        mixedsig = mixedsig' * nt;
        CovEvals = CovEvals / npix;

        percentvar = 100*sum(CovEvals)/covtrace;
        fprintf([' First ',num2str(nPCs),' PCs contain ',num2str(percentvar,3),'%% of the variance.\n'])
    end

    function [mixedfilters] = reload_moviedata(npix, mov, mixedsig, CovEvals)
        %-----------------------
        % Re-load movie data
        nPCs = size(mixedsig,1);

        Sinv = inv(diag(CovEvals.^(1/2)));

        movtm = mean(mov,1); % Average over space
        movuse = mov - ones(npix,1) * movtm;
        mixedfilters = reshape(movuse * mixedsig' * Sinv, npix, nPCs);
    end

    function j = tiff_frames(fn)
        %
        % n = tiff_frames(filename)
        %
        % Returns the number of slices in a TIFF stack.
        % MJLM: Changed to use imfinfo
        %
% 		j = numel(imfinfo(fn));
        j = sizetiff_RR(fn);
    end
end