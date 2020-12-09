function main_imaging_real_data_multiple_dataSets(param_global)


%% reading params
% paths
if ~isfield(param_global, 'pathData'),param_global.pathData = [pwd,filesep,'data',filesep];
    disp(['Assuming data are in : ',param_global.pathData ]); end
if ~isfield(param_global, 'pathResults'),param_global.pathResults = [pwd,filesep,'imaging_results',filesep];
    disp(['Results are saved in : ',param_global.pathResults ]); end

% data specific [CRUCIAL]
if ~isfield(param_global, 'visibilityFileName'),  error('No data available'); end
if ~isfield(param_global, 'ObsFreq')
    fprintf('\nCRITICAL: observation frequency is not available\n')
    param_global.ObsFreq = 4536E6;
end %in arcsec

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
nDataSets = length(visibilityFileNameShort);
visibilityFileName  = cell(nDataSets,1);
for iDataSets = 1:nDataSets
    visibilityFileName{iDataSets} = [pathData,'/', visibilityFileNameShort{iDataSets}];
end
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
nMeasPerCh = zeros(nDataSets,1);

dataCells = cell(nDataSets,1);
GCells    = cell(nDataSets,1);
precondWCells   = cell(nDataSets,1);
FourierIdxWCells = cell(nDataSets,1);
epsilonVect = cell(nDataSets,1);
for iDataSets =1 :nDataSets
    [dataVect, ucorr, vcoor,wcoor, nWw,timeVect,pixelSize] = util_load_real_data(visibilityFileName{iDataSets}, param_real_data);
    nMeasPerCh(iDataSets) = length(ucorr);
    %% FoV info
    if iDataSets ==1
       FoVx = pixelSize * Nx *pi/180/3600;
       FoVy = pixelSize * Ny *pi/180/3600;
       uvGridSizex   = 1/(nufft.ox*FoVx);
       uvGridSizey   = 1/(nufft.oy*FoVy);
       minGridSize   = min(uvGridSizey,uvGridSizex); %smallest size of the gridcell in the spatial Fourier domain
    end
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
        if ~isempty(nDataBlk),   szDataBlk = floor(nMeasPerCh(iDataSets)/nDataBlk);
        else, szDataBlk = szDataBlk*(nMeasPerCh(iDataSets)>=2*szDataBlk) + nMeasPerCh(iDataSets) * (nMeasPerCh(iDataSets)< 2*szDataBlk);
        end
        param_block.size = szDataBlk;
        param_block.pos = nMeasPerCh(iDataSets); % if single data set
        dummy_block = util_time_based_block_sp_ar(timeVect,param_block);
        dataPartitions = dummy_block.block;
        param_block_structure.use_manual_partitioning = doSnapshotsBlocking;
        param_block_structure.use_equal_partitioning  = 0;
    else
        dataPartitions =nMeasPerCh;
        param_block_structure.use_equal_partitioning = 1;
        param_block_structure.equal_partitioning_no = floor(nMeasPerCh/szDataBlk);
        param_block_structure.use_manual_partitioning=0;
    end
    param_block_structure.use_density_partitioning = 0;
    param_block_structure.use_manual_frequency_partitioning = 0;
    param_block_structure.use_uniform_partitioning = 0;
    param_block_structure.partition = dataPartitions;
    
    [u, v, ~, uvidx, precondWCells{iDataSets}, nW] = util_gen_block_structure(ucorr, vcoor, aWw, nWw, param_block_structure); % this function should include try catch statements
    %
    nBlocks = length(v);
    dataCells{iDataSets} = cell(nBlocks, 1);
    for jBlk = 1:nBlocks
        dataCells{iDataSets}{jBlk} = dataVect(uvidx{jBlk});
    end
    %clear unnecessary vars at this stage
    clear ucoor vcoor timeVect;
    %% measurement operator initialization
    fprintf('\nInitializing the NUFFT operator');
    tstart = tic;
    [A, At, GCells{iDataSets}, FourierIdxWCells{iDataSets}, ~] = op_p_nufft([v u], [Ny Nx], [nufft.Ky nufft.Kx],...
        [nufft.oy*Ny nufft.ox*Nx], [Ny/2 Nx/2], nW);
    tend = toc(tstart);
    %clear unnecessary vars at this stage
    clear u v  timeVect nW;

    fprintf('Initialization runtime: %ds\n', ceil(tend));
    
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
    for jBlk = 1:nBlocks
        GBis = sparse(size(GCells{iDataSets}{jBlk},1), nFourier);
        GBis(:,FourierIdxWCells{iDataSets}{jBlk}) = GCells{iDataSets}{jBlk};
        GCells{iDataSets}{jBlk} = GBis;
    end
    if  doWProjection
        fprintf('\nINFO: w-correction is enabled \n')
        param_wterm.FoV =[FoVy; FoVx];
        param_wterm.ox = [nufft.oy ;nufft.ox];
        param_wterm.gImDims = [Ny; Nx];
        % update each block of the G matrix
        for jBlk = 1:nBlocks
            GCells{iDataSets}{jBlk} =  getWprojGmatrix(GCells{iDataSets}{jBlk},...
                wcoor(uvidx{jBlk}),param_wterm,wproj.CEnergyL2,wproj.GEnergyL2);
        end        
    end
    %clear unnecessary vars at this stage
    clear wcoor
   
    %% Set l2 constraints
    epsilonVect{iDataSets} = zeros(nBlocks,1);
    if doGenerateEpsNNLS
        %Run NNLS for to initialise the l2 bound
        param_nnls.verbose = 0.5; % print log or not
        param_nnls.rel_obj = 1e-3; % stopping criterion
        param_nnls.max_iter = 1000; % max number of iterations
        param_nnls.sol_steps = [inf]; % saves images at the given iterations
        param_nnls.beta = 1;
        partitionBis = 1;
        for jBlk = 1 : nBlocks
            partitionBis = [partitionBis size(GCells{iDataSets}{jBlk},1)];
        end
        % solve nnls per block
        fprintf('\nEstimating L2 bounds --> solving NNLS for each data block ..\n')
        [result_nnls.sol, result_nnls.L2_v] = solver_fb_nnls(vertcat(dataCells{iDataSets}{:}), ...
            @(x) vertcat(GCells{iDataSets}{:}) * A(x),@(x) At(vertcat(GCells{iDataSets}{:})' * x),param_nnls);
        residual_nnls = vertcat(dataCells{iDataSets}{:}) - vertcat(GCells{iDataSets}{:}) * A(result_nnls.sol);
        
        for jBlk = 1 : nBlocks
             epsilonVect{iDataSets}(jBlk)  = norm(residual_nnls(sum(partitionBis(1:jBlk)):sum(partitionBis(1:jBlk+1))-1));
        end
        clear *_nnls
    else
        fprintf('\nEstimating L2 bounds --> assuming Chi square distribution\n')
        for jBlk = 1 : nBlocks
             epsilonVect{iDataSets}(jBlk) = sqrt(numel(dataCells{iDataSets}{jBlk}) +2*sqrt( numel(dataCells{iDataSets}{jBlk}) ) );
        end
    end
    
end
%% combine data and measurement operator
dataCells = vertcat(dataCells{:});
dataVect  =cell2mat(dataCells);

GCells    = vertcat(GCells{:});
precondWCells = vertcat(precondWCells{:});
FourierIdxWCells = vertcat(FourierIdxWCells{:});
%% Set l2 constraints
epsilonVect = vertcat(epsilonVect{:});

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
%% compute the spectral norm of the measurement operator
fprintf('Computing operator norms ...\n');
dummyG = cell2mat(GCells);
fprintf('Natural W ...\n');
opNorm = op_norm(@(x) dummyG * A(x), @(x) At(dummyG' * x), [Ny, Nx], 1e-6, 200, verbosity);

fprintf('Precond. W ...\n');
dummyPrecondWCells = sqrt(cell2mat(precondWCells));
opPrecondNorm = op_norm(@(x) dummyPrecondWCells .* (dummyG * A(x)), ...
    @(x) At(dummyG' * (dummyPrecondWCells .* x)), ...
    [Ny, Nx], 1e-6, 500, verbosity);
fprintf('INFO:Preconditioning enabled. Operator''s spectral norm: %f ...\n',opPrecondNorm);
%clear unnecessary vars at this stage
clear dummyG dummyPrecondWCells;

%% Noise estimate in the sparsity basis
bwOp = @(x) real(At(cell2mat(GCells)' * x));
fwOp  = @(x) cell2mat(GCells) * A(x);
nMeasPerCh = sum(nMeasPerCh) ;
dirac = zeros(Ny,Nx);
dirac(Ny/2 +1,Nx/2 +1) =1;
peakPSF = max(max(bwOp(fwOp(dirac))));
noise = (randn(nMeasPerCh,1)+1i*(randn(nMeasPerCh,1)))./sqrt(2);
noise_map = bwOp(noise);
%noiseLevelDict = std(noise_map(:)./opNorm);
dirty = bwOp(cell2mat(dataCells))./peakPSF;
noiseLevelDict = norm(noise_map(:)./peakPSF);
fprintf('\nINFO: Noise level in wavelet space %f\n ',noiseLevelDict)
clear  noise  noise_map  dirac bwOp fwOp ;
%% rearrange the G matrix
for jBlk = 1:length(GCells)
    FourierIdxWCells{jBlk} = false(nFourier, 1);
    GBis =  GCells{jBlk};
    FourierIdxWCells{jBlk} = any( GBis, 1).';
    GCells{jBlk}= GBis(:, FourierIdxWCells{jBlk});
    GBis = [];
end

%% sparsity operator initialization
[Psi, Psit] = op_p_sp_wlt_basis(wvlt.basis , wvlt.nlevel, Ny, Nx);
[Psiw, Psitw] = op_sp_wlt_basis(wvlt.basis , wvlt.nlevel, Ny, Nx);
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
param_pdfb_precond.reweight_min_steps_rel_obj = 100;
param_pdfb_precond.reweight_max_reweight_itr = param_pdfb_precond.max_iter - 100;% max(500,param_pdfb_precond.max_iter - 100);
param_pdfb_precond.reweight_alpha = 10*noiseLevelDict *ones(1,wvlt.nDictionnaries) ;
param_pdfb_precond.reweight_alpha_ff = (0.5*(sqrt(5)-1));
param_pdfb_precond.reweight_abs_of_max = inf;
param_pdfb_precond.total_reweights = 20;
param_pdfb_precond.reweightBound = noiseLevelDict;

%--adaptive Eps -%
param_pdfb_precond.use_adapt_bound_eps = doAdaptiveEpsilonUpdate;
param_pdfb_precond.adapt_bound_steps = 100;
param_pdfb_precond.adapt_bound_rel_obj = 1e-4;
param_pdfb_precond.adapt_bound_tol = 1e-2;
param_pdfb_precond.adapt_bound_start = 250;


%% run adaptive PPD to compute the solution
fprintf('Starting algorithm:\n\n');
dwtmode('per');
tstart_a = tic;
fprintf('Running Adaptive PPD \n');
[result_st.sol, result_st.L1_v, result_st.L1_vp, ...
    result_st.L2_v,result_st.L2_vp, result_st.delta_v,...
    result_st.sol_v, result_st.no_sub_itr_v, ~, ~, ...
    result_st.sol_reweight_v] ...
    = solver_adaptive_ppd_fb(dataCells,...
    epsilonT, epsilonTs, epsilon, epsilons, ...
    A, At, GCells, precondWCells, FourierIdxWCells, ...
    Psi, Psit, Psiw, Psitw, ...
    param_pdfb_precond);

tend = toc(tstart_a);
result_st.runtime = tend;
result_st.no_itr = length(result_st.L1_v);
try
    for jBlk = 1:length(GCells)
        GBis = sparse(size(GCells{jBlk},1), nFourier);
        GBis(:,FourierIdxWCells{jBlk}) = GCells{jBlk};
        GCells{jBlk} = GBis;
    end
    GMatrix = cell2mat(GCells);    
    result_st.residualImage= real(At(GMatrix'*(dataVect-GMatrix*(A(result_st.sol)))));
end

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
    try

    fitsResName = [pathResults,resultFileName,'_RESIDUAL.fits'];
    fitswrite(result_st.residualImage,fitsResName);
    fitsResName = [pathResults,resultFileName,'_RESIDUAL_NZ.fits'];
    fitswrite(result_st.residualImage./peakPSF,fitsResName);
    end
    fitsDirtyName = [pathResults,resultFileName,'_DIRTY.fits'];
    fitswrite(dirty,fitsDirtyName);
    
    save([pathResults,resultFileName,'.mat'],'result_st','param_pdfb_precond','-v7.3')
end
fprintf(' Algorithm runtime: %ds\n\n', ceil(tend));

