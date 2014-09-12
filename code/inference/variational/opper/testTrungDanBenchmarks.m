function  testTrungDanBenchmarks(  )
%TESTTRUNGDANBENCHMARKS Test Trung's code on Dan's benchmarks
%   For single output, Trung's code is an improved version of Oppers
% This code (and Trung's code) uses the new version of gpml

clear all; clc; close all;
DATADIR = 'dataDan';
benchmarks = {'lineardata'};


for i = 1 : length(benchmarks)
  evaluateBenchmark(DATADIR, benchmarks{i});

end
 


return;


%% Evaluate a single benchmark
function evaluateBenchmark(DATADIR, benchmark)
% Just avoids Matlab sending me stupid warning
f = [];  dfunc = []; noise = [];  train = [];
test = [];  x = [];  y = [];

        
load([DATADIR, '/', benchmark], 'f', 'dfunc', 'func', 'noise', 'train', ...
            'test', 'x', 'y');
train = train+1; test = test+1; % shifts indices to Matlab 
[nFolds, nTrain] = size(train);
strFunc = parsePython2Matlab(func);
eval(['fwdFunc = @(theta)',strFunc, ';']);
for k = 1 : nFolds
    [xtrain, ftrain, ytrain, xtest, ftest, ytest] = ...
        readSingleFold(x, f, y, train, test,k);
    
   [mufPred, sigmafStar, yStar, pred] = runTrungSingle(fwdFunc, xtrain, ytrain, xtest, ...
                                                        noise);
    
   evaluateLatentPredictions(xtrain, ftrain, xtest, mufPred, sigmafStar, ...
                             ftest);                                  
                         
  evaluateObservablePredictions(xtrain, xtest, ytest, yStar);   
  
  fname = ['resultsOpper/', benchmark, '_k', num2str(k), '.mat'];
  %save(fname, 'mufPred', 'sigmafStar', 'yStar', 'pred');
end



%% Runs Trung's code on a Single Problem
function [mufStar, sigmafStar, yStar, pred] = runTrungSingle(fwdFunc, xtrain, ytrain, xtest,...
                                                        noise)
% lets assume that we are given the noise for now
% We can have it as an additional parameter of the covariance
% May be not,  I think I need to learn it
[conf m] = getConfigTrung(fwdFunc);

%% Learn Posterior Parameters
 m                = learnFullGaussian(m,conf);
[pred.mu, pred.sigma] = feval(m.pred, m, conf, m.xt); 
pred.Sigma = diag(pred.sigma); % we don't really care about cross-covariances
%% I AM HERE
% We need to use my functions may be to predict yStsr
 
 param  = opperMainAll(xtrain,ytrain, param, conf);
pred   = opperPredictLatent( param.meanFunc, param.covFunc, param.loghyper, ...
                           param.prior, param.post, xtrain, xtest);                                     
mufStar    = pred.mu;
sigmafStar = sqrt(diag(pred.Sigma));

yStar     = opperPredictObservable(pred, param.fwdFunc, param.fwd, conf );
                                   
return;





function [conf m] = getConfigTrung(fwdFunc, xtrain, ytrain, xtest)

covfunc = @covMatern5iso;

%% Configuration
conf.nsamples               = 2000;      % Number of samples for MC estimates
conf.covfunc                = covfunc;
conf.maxiter                = 200;
conf.displayInterval        = 20;
conf.checkVarianceReduction = false;
conf.learnhyp               = true;


%% Model 
myLogLike = @(y,f,hyp) llhGaussian(y,feval(fwdFunc,f),hyp);
[D N] =  size(xtrain);
m.x   = xtrain; m.y = ytrain; m.xt = xtest;
m.N   = N; 
m.Q   = 1; % Only a single latent function
m.pars.M = ytrain;                % Variational
m.pars.L = log(ones(m.N*m.Q,1));  % parameters
m.pars.hyp.covfunc = covfunc;   % covariance hyperparameters
m.pars.hyp.cov = cell(m.Q,1);
m.pars.hyp.cov{1}  = log(ones(eval(feval(covfunc)),1));
m.likfunc          = myLogLike;
m.pred             = @predRegression;
m.pars.hyp.lik     = 0.5*log(var(ytrain,1)/4 + 1e-4);
m.pars.hyp.likfunc = m.likfunc;
  


return;


%% covMatern5sio
function K = covMatern5iso(varargin)

K = covMaterniso(5, varargin);

return;
  
  
