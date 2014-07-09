%testGaussianLinearGaussian.m
% Gaussian Prior, Linear Fwd model, Gaussian Likelihood
%   Detailed explanation goes here
%clc;
clear all;  
close all;
rng('default');

%method = @varLinearIterativeGaussAlternate;
method = @varLinearGaussStochastic;

%% Get data and parameters of generation process
N = 100; % Number of observations
M = 10; % Number of parameters
[y param]  = generateDataLinear(M, N);
param.post = {};


%% sets configuration
conf.maxiter = 500;
conf.tol     = 1e-3;
conf.fwdfunc = @(xx) linearFwdModel(xx, param.fwd{:});

%% These are settings for stochastic method
conf.alpha   = 0.01; % learning rate
pbatch       = 1;    % proportion of datapoints in each batch
conf.nbatch  = round(pbatch*length(y));


%% Runs variational
[ param, muPred, sigmaPred ] = feval(method, y, param, conf);


%% True posterior
[muTrue, SigmaTrue] = getExactPosteriorLinearGaussian(param,y);


%% Prints error results
figure;
scatter(muTrue, muPred); xlabel('True Mu'); ylabel('Pred Mu');
set_equal_axes();
figure;
scatter(diag(SigmaTrue), sigmaPred); xlabel('True var'); ylabel('Pred var');
set(gca, 'Xscale', 'log'); set(gca, 'Yscale', 'log');
set_equal_axes();



