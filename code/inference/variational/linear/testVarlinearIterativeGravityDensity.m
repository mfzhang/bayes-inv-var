% testVarlinearIterativeGravityDensity.m
% test iterative closed-form updates for linearization in variational
% inference
clear all; clc;

addpath(genpath('~/Dropbox/Matlab/DERIVESTsuite')); % For Jacobian computations
addpath(genpath('~/Dropbox/Matlab/export_fig')); % for exporting figures
addpath(genpath('~/Dropbox/Matlab/gpml-matlab')); % solve_chol


path(path(), '/Users/edwinbonilla/Documents/home/Research/projects/geothermal/github_repo/geotherML/inversion');
path(path(), genpath('~/Dropbox/Matlab/minFunc_2012'));
path(path(), genpath('~/Dropbox/Matlab/utils'));

%GRAVDIR = './gravdata/';
GRAVDIR = ''; % gravdata folder should be in the path
EXAMPLE = 1; % anomaly
str_suffix1    = ['example','_',num2str(EXAMPLE),'_dx2_dy2_dz1.mat'];
str_suffix    = ['example','_',num2str(EXAMPLE),'_dx2_dy2_dz1_ptrain_0.1.mat'];
data_fname     = [GRAVDIR, 'test_data_',num2str(EXAMPLE),'.txt'];
G_fname        = [GRAVDIR, 'G_synthetic_', str_suffix];
grid_fname     = [GRAVDIR, 'grid_synthetic_', str_suffix1];
results_fname  = [GRAVDIR, 'results_varlinear_example', '_', str_suffix ];
LOAD_RESULTS   = 0;
MAXITER        = 1000; % Maximum number of iterations
TOL            = 1e-6;

%% loads gravity observations
idx = []; load(G_fname, 'idx');
data = load(data_fname);
data(:,1:3)  = data(:,1:3)/1000; % m to km
xx    = data(idx,1); 
yy    = data(idx,2); 
zz    = data(idx,3); %  
dobs = data(idx,4); % gravity observations
scatter(xx, yy, [], dobs);

%% loads sensitivity matrix 
load(G_fname, 'G');

  
%% algorithm configuration
rbconf.fwdfunc        =  @(xx)simple_grav_fwd(xx,G); % forward model
rbconf.D              =   size(G,2); % Number of parameters

%% Canonical parameters of prior distribution: nu and precision
mu0         = zeros(rbconf.D,1); % prior is zero mean unit variance
sigma0      = 1;
lambda0     = 1/sigma0;                 % isotropic gaussian
Lambda0     = lambda0*eye(rbconf.D);
nu0         = Lambda0*mu0;
param.prior = {nu0; Lambda0};

%% Parameters of likelihood
sigmay = 1e-7; % assume low nois observations
lambday = 1./sigmay;
param.like = {lambday};

%% empty posterior parameters force default initialization
param.post = {};

rbconf.maxiter = MAXITER;
rbconf.tol    = TOL;
%[ param, muPred, sigmaPred ] = varLinearIterativeGauss(dobs, param, rbconf);
[ param, muPred, sigmaPred ] = varLinearIterativeGaussAlternate(dobs, param, rbconf);


%% loads variables useful for plotting results
load(grid_fname, 'L_X', 'L_Y', 'L_Z', 'u_x', 'u_y', 'u_z');

%% Visualize the predicted anomaly
PRED_RHO = reshape(muPred, L_X, L_Y, L_Z);
[X,Y,Z] = meshgrid(u_x, u_y, u_z);
X = X(1:L_Y, 1:L_X, 1:L_Z); Y = Y(1:L_Y, 1:L_X, 1:L_Z); Z = Z(1:L_Y, 1:L_X, 1:L_Z);
figure,visualise_anomaly_local(X, Y, Z, PRED_RHO, []);
title('Predicted Anomaly by GP Inversion');


%% Visualize the variance
PRED_RHO = reshape(sigmaPred, L_X, L_Y, L_Z);
[X,Y,Z] = meshgrid(u_x, u_y, u_z);
X = X(1:L_Y, 1:L_X, 1:L_Z); Y = Y(1:L_Y, 1:L_X, 1:L_Z); Z = Z(1:L_Y, 1:L_X, 1:L_Z);
figure, visualise_anomaly_local(X, Y, Z, PRED_RHO, []);
title('Variance of Predicted Anomaly by GP inversion');
% 
% 










