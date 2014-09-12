% testOpperGP.m
% Test Opper's method in the GP Case
clear all; clc; close all;

%% DELETE ME
%rng('default');


%% Generate Data from GP prior and Linear Fwd Model
N       = 100; % Number of observations
PTRAIN  = 0.1;
[data trueParam] = generateDataGPPriorLinearFwd( N, PTRAIN );


%% We use true param that generatd the data
param          = trueParam;
param.post     =  {};
param.loghyper = [];

%% configuration
conf.MaxIter       = 1000;
conf.alpha         = 0.01; % learning rate
conf.S             = 10000; % Number of samples for MC estimate
conf.elogMethod    = 'mc'; % 'quad' for quadrature or 'mc' for montecarlo

%% Runs Opper's Algorithmn
[ param, muPred, sigmaPred, nelbo ] = opperMainAll(data.xtrain, data.ytrain, ...
                                     param, conf);


%% True posterior
[muTrue, SigmaTrue] = getExactPosteriorLinearGP(trueParam, data.xtrain, data.ytrain);
% 
% 
%% Prints error results
figure;
scatter(muTrue, muPred); xlabel('True Mu'); ylabel('Pred Mu');
set_equal_axes();
figure;
scatter(diag(SigmaTrue), sigmaPred); xlabel('True var'); ylabel('Pred var');
set(gca, 'Xscale', 'log'); set(gca, 'Yscale', 'log');
set_equal_axes();


%% Predicts latent functions
pred   = opperPredictLatent( param.meanFunc, param.covFunc, param.loghyper, ...
                           param.prior, param.post, data.xtrain, data.xtest);                       
muStar    = pred.mu;
sigmaStar = sqrt(diag(pred.Sigma));
figure;
plot_confidence_interval(data.xtest,muStar,sigmaStar,1.96);
hold on;
plot(data.xtest, data.ftest, 'g--', 'MarkerSize', 12);
plot(data.xtrain, data.ftrain, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
title('Latent Function');
legend({'Pred Error Bars', 'Pred Mean', 'True Test', 'Training'});

%% predicts obervable variable
conf.elogMethod    = 'mc'; % 'quad' for quadrature or 'mc' for montecarlo
conf.S             = 10000;
yStar = opperPredictObservable(pred, param.fwdFunc, ...
                                       param.fwd, conf );
figure;
plot(data.xtest, yStar, 'k-', 'LineWidth', 2, 'MarkerSize', 12);
hold on;
plot(data.xtest, data.ytest, 'g--', 'MarkerSize', 12);
hold on;
plot(data.xtrain, min(data.ytest)*ones(size(data.xtrain,1),1), 'rx', 'MarkerSize', 12);

title('Observable Function');
legend({'Predicted', 'True', 'Training Input Location'});

 