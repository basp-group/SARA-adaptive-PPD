function main_imaging_real_data(param_global)
%% This is a script for running adaptive PPD on real data sets
% The data used here are observations of the supernova remnant 3C391 using the VLA
% For more details please refer to following VLA tutorial
% https://casaguides.nrao.edu/index.php/VLA_Continuum_Tutorial_3C391-CASA5.0.0
%calibrated data can be found here https://casa.nrao.edu/Data/EVLA/3C391/EVLA_3C391_FinalCalibratedMosaicMS.tgz

%% reading params
% paths
if ~isfield(param_global, 'pathData'),param_global.pathData = [pwd,filesep,'data',filesep];
disp(['Assuming data are in : ',param_global.pathData ]); end
if ~isfield(param_global, 'pathResults'),param_global.pathResults = [pwd,filesep,'imaging_results',filesep]; 
disp(['Results are saved in : ',param_global.pathResults ]); end

% data specific [CRUCIAL]
if ~isfield(param_global, 'visibilityFileName'),  param_global.visibilityFileName = '3c391_I00';end
if ~isfield(param_global, 'ObsFreq'),   param_global.ObsFreq = 4536E6; end %in arcsec

% image resolution & dimensions [CRUCIAL]: default 1.5x nominal resolution
if ~isfield(param_global, 'pixelSize'),   param_global.pixelSize = []; end %in arcsec
if isempty (param_global.pixelSize),  param_global.imageResolution = 'nominal';
else,  param_global.imageResolution = 'user_defined';
end
if ~isfield(param_global, 'Nx'),  param_global.Nx = 1024; end
if ~isfield(param_global, 'Ny'),  param_global.Ny = 1024; end

% data blocks
if ~isfield(param_global, 'nDataBlk'),  param_global.nDataBlk = []; end
if ~isfield(param_global, 'sizeDataBlk'),  param_global.sizeDataBlk = []; end
if isempty(param_global.sizeDataBlk) && isempty (param_global.nDataBlk)
    param_global.sizeDataBlk = 2e5;
end
% sparsity dict.
if ~isfield(param_global,'wavelet_level' ), param_global.wavelet_level = 4;  end
if ~isfield(param_global,'wavelet_basis' ), param_global.wavelet_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'};  end
%  algo related
if ~isfield(param_global, 'l1_reg') , param_global.l1_reg   = 1e-4; end
if ~isfield(param_global, 'rel_obj') , param_global.rel_obj = 1e-6; end
if ~isfield(param_global, 'max_iter') , param_global.max_iter = 1e6; end

% w projection
if ~isfield(param_global, 'CEnergyL2') , param_global.CEnergyL2 =1-5e-5; end
if ~isfield(param_global, 'GEnergyL2') , param_global.GEnergyL2 =1-1e-4; end

% flags

if ~isfield(param_global,'flag_adaptive_l2bounds' ), param_global.flag_adaptive_l2bounds =0 ; end
if ~isfield(param_global,'flag_generate_eps_nnls' ),  param_global.flag_generate_eps_nnls = 0 ; end
if param_global.flag_adaptive_l2bounds; param_global.flag_generate_eps_nnls =1;end

if ~isfield(param_global,'flag_wProjection' ), param_global.flag_wProjection = 0 ; end
if ~isfield(param_global,'flag_compute_Anorm' ),    param_global.flag_compute_Anorm =0 ; end
if ~isfield(param_global,'flag_save_results' ), param_global.flag_save_results = 1;  end
if ~isfield(param_global,'flag_verbose' ), param_global.flag_verbose = 0;  end
if ~isfield(param_global,'flag_preconditionning' ), param_global.flag_preconditionning = 1;  end
if ~isfield(param_global,'flag_manual_partitioning' ), param_global.flag_manual_partitioning = 1;  end
if ~isfield(param_global,'flag_is_weightsFixed2sigma' ), param_global.flag_is_weightsFixed2sigma = 0;  end

%% paths, data specs and imaging params !! CRUCIAL INPUT FROM USER
pathData = param_global.pathData;
pathResults = param_global.pathResults;

% data specs
visibilityFileNameShort = param_global.visibilityFileName;% data in .mat format
visibilityFileName = [pathData,'/', visibilityFileNameShort];
% Imaging params
Nx = param_global.Nx;
Ny = param_global.Ny;
pixelSize = param_global.pixelSize;% 'in arcsec'
obsFreq   = param_global.ObsFreq;% freq in 'Hz'
imageResolutionStr = param_global.imageResolution;

fprintf('\nINFO: Observation freq: %f MHz.\n',obsFreq/1E6);
switch imageResolutionStr
    case 'nominal'
        fprintf('\nWARNING: No pixelsize provided by user --> adopting 1.5x instrumental resolution.\n')
    otherwise
        fprintf('\nINFO: Pixelsize provided by user (in arcsec by default): %f asec.\n',pixelSize);
end

isWeightsFixed2Sigma =param_global.flag_is_weightsFixed2sigma ;
if isWeightsFixed2Sigma
    fprintf('\nWARNING: reading weights from data --> SIGMA col loaded\n') ;
else
    fprintf('\nINFO: reading weights from data --> WEIGHTS col loaded\n') ;
end

%% flags
doAdaptiveEpsilonUpdate = param_global.flag_adaptive_l2bounds;
doGenerateEpsNNLS= param_global.flag_generate_eps_nnls;
doComputeAnorm   = param_global.flag_compute_Anorm;
doSaveResults = param_global.flag_save_results;
verbosity     =  param_global.flag_verbose;
doPreconditionning  = param_global.flag_preconditionning;
doSnapshotsBlocking  = param_global.flag_manual_partitioning;
doWProjection = param_global.flag_wProjection;
%% config params

% sparsity prior
wvlt.nlevel = param_global.wavelet_level ; % wavelet level
wvlt.basis  = param_global.wavelet_basis; % wavelet basis to be used, always put self in last position if used
wvlt.nDictionnaries = length(wvlt.basis);

% NUFFT parameters
nufft.ox = 2 ;% oversampling factors for nufft
nufft.oy = 2 ;% oversampling factors for nufft
nufft.Kx = 8 ;% number of neighbours for nufft
nufft.Ky = 8 ;% number of neighbours for nufft
nFourier  = prod([nufft.oy*Ny nufft.ox*Nx]) ;

% algo related
l1_reg = param_global.l1_reg;
rel_obj  = param_global.rel_obj; % stopping criterion
max_iter = param_global.max_iter ; % max number of iterations

%data splitting parameters
nDataBlk= param_global.nDataBlk ;
szDataBlk= param_global.sizeDataBlk ;

% meas. op.: w projection
wproj.CEnergyL2 = param_global.CEnergyL2 ;
wproj.GEnergyL2 = param_global.GEnergyL2 ;


%% Get data, uvw coord, weights and time

param_real_data.pixelSize = pixelSize;% 'in arcsec'
param_real_data.obsFreq   = obsFreq;% freq in 'Hz'
param_real_data.imageResolutionStr =imageResolutionStr;
param_real_data.isWeightsFixed2Sigma = isWeightsFixed2Sigma;

[dataVect, ucorr, vcoor,wcoor, nWw,timeVect,pixelSize] = util_load_real_data(visibilityFileName, param_real_data);
nMeasPerCh = length(ucorr);
%% FoV info

FoVx = pixelSize * Nx *pi/180/3600;
FoVy = pixelSize * Ny *pi/180/3600;
uvGridSizex   = 1/(nufft.ox*FoVx);
uvGridSizey   = 1/(nufft.oy*FoVy);
minGridSize   = min(uvGridSizey,uvGridSizex); %smallest size of the gridcell in the spatial Fourier domain

%% Preconditioning: compute weights for Projection into the ellipsoids

if doPreconditionning
    param_precond.N = nFourier; % number of pixels in the image
    param_precond.Nox = nufft.oy*Ny; % number of pixels in the image
    param_precond.Noy = nufft.oy*Ny; % number of pixels in the image
    param_precond.gen_uniform_weight_matrix = 1; %set weighting type
    param_precond.uniform_weight_sub_pixels = 1;
    aWw = util_gen_preconditioning_matrix(ucorr, vcoor, param_precond);
else, aWw =ones(size(ucorr));
end

%% set the blocks structure
% check time vector is read along with the data
if (numel(unique(timeVect))==1) && (doSnapshotsBlocking)
    fprintf('\nWarning: time vector not found --> default blocking strategy: equal blk size \n')
    doSnapshotsBlocking = 0;
end

if doSnapshotsBlocking
    if ~isempty(nDataBlk),   szDataBlk = floor(nMeasPerCh/nDataBlk);
    else, szDataBlk = szDataBlk*(nMeasPerCh>=2*szDataBlk) + nMeasPerCh * (nMeasPerCh< 2*szDataBlk);
    end
    param_block.size = szDataBlk;
    param_block.pos = nMeasPerCh; % if single data set
    dummy_block = util_time_based_block_sp_ar(timeVect,param_block);
    dataPartitions = dummy_block.block;
    param_block_structure.use_manual_partitioning = doSnapshotsBlocking;
    param_block_structure.use_equal_partitioning  = 0;
else
    param_block_structure.use_equal_partitioning = 1;
    param_block_structure.equal_partitioning_no = floor(MeasPerCh/szDataBlk);
    param_block_structure.use_manual_partitioning=0;
end
param_block_structure.use_density_partitioning = 0;
param_block_structure.use_manual_frequency_partitioning = 0;
param_block_structure.use_uniform_partitioning = 0;
param_block_structure.partition = dataPartitions;

[u, v, ~, uvidx, aW, nW] = util_gen_block_structure(ucorr, vcoor, aWw, nWw, param_block_structure); % this function should include try catch statements
%
nBlocks = length(v);
dataCells = cell(nBlocks, 1);
for jBlk = 1:nBlocks
    dataCells{jBlk} = dataVect(uvidx{jBlk});
end
%clear unnecessary vars at this stage
clear ucoor vcoor;
%% measurement operator initialization
fprintf('\nInitializing the NUFFT operator\n');
tstart = tic;
[A, At, GCells, W, GMatrix] = op_p_nufft_up([v u], [Ny Nx], [nufft.Ky nufft.Kx],...
    [nufft.oy*Ny nufft.ox*Nx], [Ny/2 Nx/2], nW);
tend = toc(tstart);

fprintf('Initialization runtime: %ds\n\n', ceil(tend));

%% w-correction
%:AD:  check if w correction is necessary
if ~doWProjection
    effBandwidthWterm = max(abs(max(FoVy, FoVx).*wcoor(:)));
    if effBandwidthWterm > 3*minGridSize
        doWProjection = 1;
    else
        fprintf('\nINFO: no w-orrection\n')
        doWProjection = 0;
    end
end

if  doWProjection    
    fprintf('\nINFO: w-correction is enabled \n')
    param_wterm.FoV =[FoVy; FoVx];
    param_wterm.ox = [nufft.oy ;nufft.ox];
    param_wterm.gImDims = [Ny; Nx];
    % update each block of the G matrix
    for jBlk = 1:nBlocks
        GCells{jBlk} =  getWprojGmatrix(GCells{jBlk},  wcoor(uvidx{jBlk}),param_wterm,wproj.CEnergyL2,wproj.GEnergyL2);
    end
end

%% compute the spectral norm of the measurement operator
fprintf('Computing operator norms ...\n');
if doComputeAnorm
    fprintf('Natural W ...\n');
    opNorm = op_norm(@(x) GMatrix * A(x), @(x) At(GMatrix' * x), [Ny, Nx], 1e-6, 200, verbosity);
end

opPrecondNorm = op_norm(@(x) sqrt(cell2mat(aW)) .* (GMatrix * A(x)), @(x) At(GMatrix' * (sqrt(cell2mat(aW)) .* x)), [Ny, Nx], 1e-6, 200, verbosity);
fprintf('INFO:Preconditioning enabled. Operator''s spectral norm: %f ...\n',opPrecondNorm);
%clear unnecessary vars at this stage
clear ucoor vcoor u v  timeVect nW;

%% Set l2 constraints
epsilonVect = zeros(nBlocks,1);
if doGenerateEpsNNLS
    %Run NNLS for to initialise the l2 bound
    param_nnls.verbose = 0.5; % print log or not
    param_nnls.rel_obj = 1e-4; % stopping criterion
    param_nnls.max_iter = 1000; % max number of iterations
    param_nnls.sol_steps = [inf]; % saves images at the given iterations
    param_nnls.beta = 1;
    
    % solve nnls per block
    fprintf('\nEstimating L2 bounds --> solving NNLS for each data block ..\n')
    for jBlk = 1 : nBlocks
        [result_nnls.sol{jBlk}, result_nnls.L2_v{jBlk}] = solver_fb_nnls(dataCells{jBlk}, @(x) GCells{jBlk} * A(x),@(x) At((GCells{jBlk})' * x),param_nnls);
        epsilonVect(jBlk)  = result_nnls.L2_v{jBlk}(end);
    end
else
    fprintf('\nEstimating L2 bounds --> assuming Chi square distribution\n')
    for jBlk = 1 : nBlocks
        epsilonVect(jBlk) =sqrt(numel(dataCells{jBlk}) +2*sqrt( numel(dataCells{jBlk}) ) );
    end
end
l2_ball_definition = 'value';
stopping_criterion = 'l2-ball-percentage';
param_l2_ball.stop_eps_v = epsilonVect;
param_l2_ball.val_eps_v =  param_l2_ball.stop_eps_v;
param_l2_ball.l2_ball_percentage_stop = 1.01;
use_same_stop_criterion = 1;
sigma_noise = sqrt(2);
[epsilonT, epsilonTs, epsilon, epsilons] = util_gen_L2_bounds(dataCells, ... % this function should include try catch statements
    [], sigma_noise, l2_ball_definition, stopping_criterion,...
    use_same_stop_criterion, param_l2_ball);

%% sparsity operator initialization
[Psi, Psit] = op_p_sp_wlt_basis(wvlt.basis , wvlt.nlevel, Ny, Nx);
[Psiw, Psitw] = op_sp_wlt_basis(wvlt.basis , wvlt.nlevel, Ny, Nx);

%% Noise estimate in the sparsity basis
bwOp = @(x) real(At(GMatrix' * x));
fwOp  = @(x) GMatrix * A(x);
dirac = zeros(Ny,Nx);
dirac(Ny/2 +1,Nx/2 +1) =1;
peakPSF = max(max(bwOp(fwOp(dirac))));
noise = (randn(nMeasPerCh,1)+1i*(randn(nMeasPerCh,1)))./sqrt(2);
noise_map = bwOp(noise);
noiseLevelDict = std(noise_map(:)./peakPSF);
dirty = bwOp(cell2mat(dataCells))./peakPSF;

fprintf('\nINFO: Noise level in wavelet space %f\n ',noiseLevelDict)
clear GMatrix noise  noise_map  dirac bwOp fwOp ;
%% PDFB parameter structure sent to the algorithm
%----default----%
param_pdfb_precond.im = 0; % original image, used to compute the SNR
param_pdfb_precond.verbose = 1; % print log or not
param_pdfb_precond.nu1 = 1; % bound on the norm of the operator Psi
param_pdfb_precond.nu2 = opPrecondNorm; % bound on the norm of the operator A*G
param_pdfb_precond.tau = 0.49; % forward descent step size
param_pdfb_precond.lambda0 = 1; % relaxation step for primal update
param_pdfb_precond.lambda1 = 1; % relaxation step for L1 dual update
param_pdfb_precond.lambda2 = 1; % relaxation step for L2 dual update
param_pdfb_precond.sol_steps =inf; % saves images at the given iterations
param_pdfb_precond.use_proj_elipse_fb = 1;
param_pdfb_precond.elipse_proj_max_iter = 20;
param_pdfb_precond.elipse_proj_min_iter = 1;
param_pdfb_precond.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid
%---------------%
%-user defined-%
param_pdfb_precond.rel_obj = rel_obj; % stopping criterion
param_pdfb_precond.max_iter = max_iter; % max number of iterations
param_pdfb_precond.gamma = l1_reg; % convergence parameter L1 (soft thresholding parameter)
%--re-weighting-%
param_pdfb_precond.use_reweight_eps =1;

param_pdfb_precond.reweight_rel_obj = 5e-5; % criterion for performing reweighting
param_pdfb_precond.reweight_min_steps_rel_obj = 250;
param_pdfb_precond.reweight_max_reweight_itr = max(500,param_pdfb_precond.max_iter - 500);
param_pdfb_precond.reweight_alpha = ones(1,wvlt.nDictionnaries) ;
param_pdfb_precond.reweight_alpha_ff = (0.5*(sqrt(5)-1)) * ones(1,wvlt.nDictionnaries);
param_pdfb_precond.reweight_abs_of_max = inf;
param_pdfb_precond.total_reweights = 20;
param_pdfb_precond.reweightBound = noiseLevelDict;

%--adaptive Eps -%
param_pdfb_precond.use_adapt_bound_eps = doAdaptiveEpsilonUpdate;
param_pdfb_precond.adapt_bound_steps = 100;
param_pdfb_precond.adapt_bound_rel_obj = 1e-4;
param_pdfb_precond.adapt_bound_tol = 1e-3;
param_pdfb_precond.adapt_bound_start = 200;


%% run adaptive PPD to compute the solution
fprintf('Starting algorithm:\n\n');
dwtmode('per');
tstart_a = tic;
fprintf(' Running Adaptive PPD \n');
[result_st.sol, result_st.L1_v, result_st.L1_vp, ...
    result_st.L2_v,result_st.L2_vp, result_st.delta_v,...
    result_st.sol_v, result_st.no_sub_itr_v, ~, ~, ...
    result_st.sol_reweight_v] ...
    = solver_adaptive_ppd_fb(dataCells,...
    epsilonT, epsilonTs, epsilon, epsilons, ...
    A, At, GCells, aW, W, ...
    Psi, Psit, Psiw, Psitw, ...
    param_pdfb_precond);

tend = toc(tstart_a);
result_st.runtime = tend;
result_st.no_itr = length(result_st.L1_v);
try GMatrix = cell2mat(GCells);
catch,    GMatrix = vertcat(GCells{:});
end
result_st.residualImage= real(At(GMatrix'*(dataVect-GMatrix*(A(result_st.sol)))));

fprintf(' Adaptive PPD runtime: %ds\n\n', ceil(tend));

clear GCells GMatrix;
%% save results
if doSaveResults
    resultFileName =(['APPD-rslts.',...
        '.Nx',num2str(Nx),'.Ny',num2str(Ny),...
        '.dl',num2str(param_real_data.pixelSize),...
        '.eps',num2str(result_st.L2_v(end)),...
        '.g',num2str(param_pdfb_precond.gamma),...
        '.Nitr',num2str(length(result_st.L2_v))]);
    
    fitsSolName = [pathResults,resultFileName,'_SOL.fits'];
    fitswrite(result_st.sol,fitsSolName);
    fitsResName = [pathResults,resultFileName,'_RESIDUAL.fits'];
    fitswrite(result_st.residualImage,fitsResName);
    fitsResName = [pathResults,resultFileName,'_RESIDUAL_NZ.fits'];
    fitswrite(result_st.residualImage./peakPSF,fitsResName);
    fitsDirtyName = [pathResults,resultFileName,'_DIRTY.fits'];
    fitswrite(dirty,fitsDirtyName);
    
    save([pathResults,resultFileName,'.mat'],'result_st','param_pdfb_precond','-v7.3')
end
fprintf(' Algorithm runtime: %ds\n\n', ceil(tend));

