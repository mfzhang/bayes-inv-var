% testOpper.m
%
clear all; clc;

%% DELETE ME
% rng('default');


%% Generate  data from linear model
% param are the natural parameters
N          = 10; % Number of observations
[y param]  = generateDataisoLinear(N);
% load('bad_example_data.mat', 'y', 'param'); 
param.post = {};


%% configuration
conf.MaxIter       = 1000;
conf.alpha         = 0.01; % learning rate
conf.S             = 1000; % Number of samples for MC estimate
conf.elogMethod    = 'quad'; % 'quad' for quadrature or 'mc' for montecarlo

%% Runs Opper Algorith
[ param, muPred, sigmaPred ] = opperMainGaussFullReparam( y, param, conf );
 
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

