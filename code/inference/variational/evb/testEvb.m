% test EVB (Enhanced Variational Bayes)
clear all; clc;

%% Generate  data from linear model
% param are the natural parameters
N = 100; % Number of observations
M = 10;  % Number of parameters
[y param]  = generateDataLinear(M, N);
param.post = {};


%% configuration
conf.MaxIter = 500;
conf.alpha   = 0.01; % learning rate
conf.fwdfunc = @(xx) linearFwdModel(xx, param.fwd{:});
conf.S       = 500; % Number of samples for MC estimate

%% Runs EVB
figure;
[ param, muq, sigmaPred ] = evbMainGaussFull( y, param, conf );
