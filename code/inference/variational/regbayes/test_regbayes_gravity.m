% test_regbayes_gravity.m
% test regbayes on gravity problem
% Assume Known A=G 
% learn variational parameters of factorized posterior
clear all; clc; close all;

path(path(), '/Users/edwinbonilla/Documents/home/Research/projects/geothermal/github_repo/geotherML/inversion');
GRAVDIR = '';
EXAMPLE = 1; % anomaly
str_example    = num2str(EXAMPLE);
data_fname     = [GRAVDIR, 'test_data_',str_example,'.txt'];
G_fname        = [GRAVDIR, 'G_synthetic_example_', str_example, '.mat'];
grid_fname     = [GRAVDIR, 'grid_synthetic_example_',str_example, '.mat'];
results_fname  = [GRAVDIR, 'results_example_', str_example, '.mat' ];
LOAD_RESULTS   = 1;

%% loads gravity observations
data = load(data_fname);
data(:,1:3)  = data(:,1:3)/1000; % m to km
xx    = data(:,1); 
yy    = data(:,2); 
zz    = data(:,3); %  
dobs = data(:,4); % gravity observations
scatter(xx, yy, [], dobs);

%% loads sensitivity matrix 
load(G_fname, 'G');


%% algorithm configuration
rbconf.klfunc         =  @kl_gauss_iso_diag;
rbconf.elikefunc      =  @eloglike_gauss_iso_diag;
rbconf.fwdfunc        = []; % forward model
rbconf.regfunc        = []; % Regularizer
rbconf.param2vecfunc  =  @param2vec_gauss_diag;
rbconf.vec2paramfunc  =  @vec2param_gauss_diag;
rbconf.initparamfunc  =  @initparam_gauss_diag;
rbconf.optfunc        =  @minFunc; % {'minFunc', 'fminunc'}
rbconf.D              =   size(G,2); % Number of parameters


%% optimizer configuration
optconf = get_default_optconf();

%% parameters of distributions
mu0         = zeros(rbconf.D,1); % prior is zero mean unit variance
sigma0      = 1;                 % isotropic gaussian
param.prior = {mu0, sigma0};
%
sigmay = 1e-7; % assume low nois observations
param.like  =  {G, sigmay};
param.fwd   = [];

if (LOAD_RESULTS)
    load(results_fname, 'mu_pred', 'sigma_pred');
else
    %% runs variational inference
    [post_param, nlog] = learn_regbayes(dobs, rbconf, optconf, param);
    mu_pred            = post_param{1};
    sigma_pred         = post_param{2};
    save(results_fname, 'mu_pred', 'sigma_pred');
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

