%% This is a script for running adaptive PPD on real data sets
% The data used here are observations of the supernova remnant 3C391 using the VLA 
% For more details please refer to following VLA tutorial
% https://casaguides.nrao.edu/index.php/VLA_Continuum_Tutorial_3C391-CASA5.0.0
%calibrated data can be found here https://casa.nrao.edu/Data/EVLA/3C391/EVLA_3C391_FinalCalibratedMosaicMS.tgz
%% 
clear; clc ; close all

try
        % set NUFFT path
        irt= ['./irt/setup.m'];
        project  = './';
        run(irt)
end

ppd_repository = [project,'ppd-for-ri-master/'];
addpath(ppd_repository);
addpath([ppd_repository, 'lib/']);
addpath([ppd_repository, 'alg/']);
addpath([project,'lib_appd'])

save_path = [project,'imaging_results/'];
data_path = [project,'data/'];
mkdir(save_path)
%%
% general config parameters
verbosity = 1;
use_same_stop_criterion = 1;
flag_single_data_set = 1;
flag_save_results =1;

% preconditioning
param_precond.gen_uniform_weight_matrix = 1; %set weighting type
param_precond.uniform_weight_sub_pixels = 1;

% sparsity prior
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; % wavelet basis to be used
dictn = 9;
nlevel = 4; % wavelet level

% NUFFT parameters
ox = 2 ;% oversampling factors for nufft
oy = 2 ;% oversampling factors for nufft
Kx = 8 ;% number of neighbours for nufft
Ky = 8 ;% number of neighbours for nufft

% get input data
visibility_file_name_short = '3c391_I00';% data in .mat format
visibility_file_name = [data_path,visibility_file_name_short];

%data splitting parameters
param_block.size = 20000;
param_block.snapshot=0;
param_block.pos = 68699;%  if a single data set set it to the total number or measurements else set it to [M1 M2]; M1, M2 are sizes of the concatenated data sets
param_block_structure.use_manual_partitioning = 1;
use_gridded_data = 0; % flag setting for generating gridded data

% Imaging parameteres
param_real_data.image_size_Nx = 512;
param_real_data.image_size_Ny = 512;
param_real_data.pixel_size = 2.5; % 'in arcsec'
param_real_data.freq = 4536E6;% freq in 'Hz'

%% Get data and build operators 

script_get_real_data_build_operators


%% Run NNLS for to initialise the l2 bound

param_nnls.im = zeros(param_real_data.image_size_Nx,param_real_data.image_size_Ny); % original image, used to compute the SNR
param_nnls.verbose = verbosity; % print log or not
param_nnls.rel_obj = 5e-5; % stopping criterion
param_nnls.max_iter = 500; % max number of iterations
param_nnls.sol_steps = [inf]; % saves images at the given iterations
param_nnls.beta = 1;
epsVect = zeros(1,length(blockStruct.blockNumber));
for i =1:(blockStruct.blockNumber)
[result_nnls.sol{i}, result_nnls.L2_v{i}] = solver_fb_nnls(y{i}, @(x) T{i} * A(x),@(x) At((T{i})' * x),param_nnls);
epsVect(i)  = result_nnls.L2_v{i}(end);
end
%% Initialize l2 bounds

l2_ball_definition = 'value';
stopping_criterion = 'l2-ball-percentage';

param_l2_ball.stop_eps_v = epsVect;
param_l2_ball.val_eps_v = 1.* param_l2_ball.stop_eps_v;
param_l2_ball.l2_ball_percentage_stop = 1.001;

[epsilonT, epsilonTs, epsilon, epsilons] = util_gen_L2_bounds(y, ...
                [], sigma_noise, l2_ball_definition, stopping_criterion,...
                use_same_stop_criterion, param_l2_ball);
 
%% PDFB parameter structure sent to the algorithm


param_pdfb_precond.im = 0; % original image, used to compute the SNR
param_pdfb_precond.verbose = 1; % print log or not
param_pdfb_precond.nu1 = 1; % bound on the norm of the operator Psi
param_pdfb_precond.nu2 = evl_precond; % bound on the norm of the operator A*G
param_pdfb_precond.tau = 0.49; % forward descent step size

param_pdfb_precond.lambda0 = 1; % relaxation step for primal update
param_pdfb_precond.lambda1 = 1; % relaxation step for L1 dual update
param_pdfb_precond.lambda2 = 1; % relaxation step for L2 dual update
param_pdfb_precond.sol_steps =inf; % saves images at the given iterations

param_pdfb_precond.use_proj_elipse_fb = 1;
param_pdfb_precond.elipse_proj_max_iter = 20;
param_pdfb_precond.elipse_proj_min_iter = 1;
param_pdfb_precond.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid

param_pdfb_precond.path='./itrs/';mkdir(param_pdfb_precond.path);
                                                                                   
param_pdfb_precond.rel_obj = 5e-6; % stopping criterion
param_pdfb_precond.max_iter = 10000; % max number of iterations
param_pdfb_precond.gamma = 0.8e-4; % convergence parameter L1 (soft thresholding parameter)

param_pdfb_precond.use_reweight_eps = 1;

reweight = param_pdfb_precond.use_reweight_eps;
param_pdfb_precond.reweight_rel_obj = 5e-5; % criterion for performing reweighting
param_pdfb_precond.reweight_min_steps_rel_obj = 250;
param_pdfb_precond.reweight_max_reweight_itr = max(500,param_pdfb_precond.max_iter - 500);
param_pdfb_precond.reweight_alpha = ones(1,dictn) ;
param_pdfb_precond.reweight_alpha_ff = 0.8 * ones(1,dictn);
param_pdfb_precond.reweight_abs_of_max = inf;
param_pdfb_precond.total_reweights = 10;
                                                                                   

param_pdfb_precond.use_adapt_bound_eps = 1;
param_pdfb_precond.adapt_bound_steps = 100;
param_pdfb_precond.adapt_bound_rel_obj = 1e-4;
param_pdfb_precond.adapt_bound_tol = 1e-3;
param_pdfb_precond.adapt_bound_start = 200;


%% run adaptive PPD to compute the solution
fprintf('Starting algorithm:\n\n');
tstart = tic;
script_run_adaptive_ppd;
tend = toc(tstart);

%% save results
if flag_save_results
result_file_name =(['AdaptivePPD-rslts.nRW.',...
num2str(param_pdfb_precond.reweight_abs_of_max),...
'.wavelet',num2str(nlevel),...
'.Nx',num2str(param_real_data.image_size_Nx),...
'.dl',num2str(param_real_data.pixel_size),...
'.eps',num2str(result_st.L2_v(end)),...
'.g',num2str(param_pdfb_precond.gamma),...
'.Nitr',num2str(length(result_st.L2_v))]);

fitsSolName = [save_path,result_file_name,'_SOL.fits'];
fitswrite(result_st.sol,fitsSolName);
fitsResName = [save_path,result_file_name,'_RESIDUAL.fits'];
fitswrite(result_st.residualImage,fitsResName);

save([save_path,result_file_name,'.mat'],'result_st','param_pdfb_precond','-v7.3')
end
fprintf(' Algorithm runtime: %ds\n\n', ceil(tend));

