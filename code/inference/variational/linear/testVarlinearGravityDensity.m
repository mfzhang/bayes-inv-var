%test_variational_linear.m
% tests linearization in variational inference

clear all; clc; close all;
addpath(genpath('~/Dropbox/Matlab/DERIVESTsuite')); % For Jacobian computations
addpath(genpath('~/Dropbox/Matlab/export_fig')); % for exporting figures


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
rbconf.klfunc         =  @kl_gauss_iso_diag;
rbconf.elikefunc      =  @eloglike_linear_gauss_isodiag;
rbconf.fwdfunc        =  @simple_grav_fwd; % forward model
rbconf.param2vecfunc  =  @param2vec_gauss_diag;
rbconf.vec2paramfunc  =  @vec2param_gauss_diag;
rbconf.initparamfunc  =  @initparam_gauss_diag;
rbconf.optfunc        =  @minFunc; % {'minFunc', 'fminunc'}
rbconf.D              =   size(G,2); % Number of parameters


%% optimizer configuration
optconf = get_optconf_varlinear();

%% parameters of prior distribution
mu0         = zeros(rbconf.D,1); % prior is zero mean unit variance
sigma0      = 1;                 % isotropic gaussian
param.prior = {mu0, sigma0};

%% Parameters of likelihood
sigmay = 1e-7; % assume low nois observations
param.like  =  {sigmay};

%% Parameters of fwd model
param.fwd   = {G};

%% empty posterior parameters force default initialization
param.post = {};


if (LOAD_RESULTS)
    load(results_fname, 'mu_pred', 'sigma_pred');
else
 
        %% runs variational inference
        [post_param, nlog] = learn_varlinear(dobs, rbconf, optconf, param);
        mu_pred            = post_param{1};
        sigma_pred         = post_param{2};
end




%% loads variables useful for plotting results
load(grid_fname, 'L_X', 'L_Y', 'L_Z', 'u_x', 'u_y', 'u_z');


%% Visualize the predicted anomaly
PRED_RHO = reshape(mu_pred, L_X, L_Y, L_Z);
[X,Y,Z] = meshgrid(u_x, u_y, u_z);
X = X(1:L_Y, 1:L_X, 1:L_Z); Y = Y(1:L_Y, 1:L_X, 1:L_Z); Z = Z(1:L_Y, 1:L_X, 1:L_Z);
figure,visualise_anomaly_local(X, Y, Z, PRED_RHO, []);
title('Predicted Anomaly by GP Inversion');

%% Visualize the variance
PRED_RHO = reshape(sigma_pred, L_X, L_Y, L_Z);
[X,Y,Z] = meshgrid(u_x, u_y, u_z);
X = X(1:L_Y, 1:L_X, 1:L_Z); Y = Y(1:L_Y, 1:L_X, 1:L_Z); Z = Z(1:L_Y, 1:L_X, 1:L_Z);
figure, visualise_anomaly_local(X, Y, Z, PRED_RHO, []);
title('Variance of Predicted Anomaly by GP inversion');

s = datestr(now);
save(s);



