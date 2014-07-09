% test_regbayes_gravity_reg.m
% test regbayes on gravity problem with squared regularizer term
% it also learns the sensitivity of the fwd model
clear all; clc; close all;

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
results_fname  = [GRAVDIR, 'results_example', '_', str_suffix ];
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
rbconf.elikefunc      =  @eloglike_gauss_iso_diag;
rbconf.fwdfunc        =  @simple_grav_fwd; % forward model
rbconf.regfunc        =  @sqreg_gauss_diag; % Regularizer
rbconf.param2vecfunc  =  @param2vec_gauss_diag;
rbconf.vec2paramfunc  =  @vec2param_gauss_diag;
rbconf.initparamfunc  =  @initparam_gauss_diag;
rbconf.optfunc        =  @minFunc; % {'minFunc', 'fminunc'}
rbconf.D              =   size(G,2); % Number of parameters


%% optimizer configuration
optconf = get_default_optconf();

%% parameters of prior distribution
mu0         = zeros(rbconf.D,1); % prior is zero mean unit variance
sigma0      = 1;                 % isotropic gaussian
param.prior = {mu0, sigma0};

%% Parameters of likelihood
sigmay = 1e-7; % assume low nois observations
A = rand(length(dobs), rbconf.D);  % sensitivity matrix of approximate likelihood
% TODO: we can initialize A above by fitting a few samples from fwd model
param.like  =  {A, sigmay};

%% Parameters of fwd model
param.fwd   = {G};

%% paramters of regularizer
lambda_reg  = 1;
param.reg   = {lambda_reg};

%% empty posterior parameters force default initialization
param.post = {};

if (LOAD_RESULTS)
    load(results_fname, 'mu_pred', 'sigma_pred');
else
    %% Here we alternate between learning (for likelihood) and learning
    % variational parameters
    nlog = zeros(100,1);
    for i = 1 : 1000
    
        %% runs variational inference
        [post_param, nlog(i)] = learn_regbayes(dobs, rbconf, optconf, param);
        mu_pred            = post_param{1};
        sigma_pred         = post_param{2};
        %save(results_fname, 'mu_pred', 'sigma_pred');
        
       % showing nelbo value here
       % theta = feval(rbconf.param2vecfunc, post_param);
       % nlog  =  regbayes_negelbo(theta, dobs, rbconf, param) 
        
        param.like{1} = update_sensitivity_gauss(dobs, param.like{2}, lambda_reg, ...
                        rbconf.fwdfunc, param.fwd, mu_pred, diag(sigma_pred));

       %nlog  =  regbayes_negelbo(theta, dobs, rbconf, param) 
       %pause;
       
       
       %% Updates posterior parameteres
       param.post{1} = mu_pred;
       param.post{2} = sigma_pred;
    end
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

save(['regbayes-results-', s]);

