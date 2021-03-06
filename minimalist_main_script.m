%% This is a script for running adaptive PPD on real data sets
% The data used here are observations of the supernova remnant 3C391 using the VLA
% For more details please refer to following VLA tutorial
% https://casaguides.nrao.edu/index.php/VLA_Continuum_Tutorial_3C391-CASA5.0.0
%calibrated data can be found here https://casa.nrao.edu/Data/EVLA/3C391/EVLA_3C391_FinalCalibratedMosaicMS.tgz
clear ; clc; close all;
%% setting up paths
pathProject  = pwd;
if pathProject(end) ~= filesep % make sure there is a '/' at end of directory
    pathProject = [pathProject filesep];
end

pathLibNUFFT = [pathProject,'lib',filesep, 'measurement-operator',filesep,'nufft',filesep];
addpath(pathLibNUFFT);
% Adaptive PPD: Dabbech et al. 2018
pathLibAPPD=[pathProject,'src',filesep];
addpath(genpath(pathLibAPPD));

% PPD: Onose et al 2016
pathLibPPD = [pathProject,'lib',filesep,'measurement-operator',filesep, 'utils',filesep];
addpath(genpath(pathLibPPD));

% W-correction: dabbech et al 2017
pathWCorrection = [pathProject,'lib',filesep,'measurement-operator',filesep, 'lib', filesep,'wproj_utilities',filesep];
addpath(pathWCorrection);

% DATA dir.
pathData = [pathProject,'data',filesep];

% Results dir.
pathResults = [pathProject,'imaging_results',filesep];
try mkdir(pathResults) 
end
%% param struct
param_global.pathData = pathData; %data path
param_global.pathResults = pathResults; %results path 

param_global.Nx = 512; %image dim.
param_global.Ny = 512; %image dim.

param_global.pixelSize = 1.98; % 'in arcsec'
param_global.ObsFreq = 1e9;% freq in 'Hz'

param_global.visibilityFileName = {'data_uv_arwa_3c353_6'} ;%supports multiple datasets, e.g.  {'data1.mat','data2.mat'}
param_global.flag_is_weightsFixed2sigma = 1; % 1: 'weights' in the data file corresponds to SIGMA. 0: 'weights' correspond to the inverse of the VARIANCE

%solver
param_global.l1_reg = 1e-4;% soft-thresholding param. of SARA
param_global.max_iter = 10000;
param_global.rel_obj = 5e-6; 

param_global.flag_adaptive_l2bounds =0 ; % no adaptive epsilon if noise statistics are known (always the case for sims)
param_global.sizeDataBlk = 200000; % not compulsory: set data block sizes in the solver


%% run Adaptive PPD
main_imaging_real_data_multiple_dataSets(param_global);



