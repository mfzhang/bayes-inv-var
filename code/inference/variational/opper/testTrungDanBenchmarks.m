function  testTrungDanBenchmarks(  )
%TESTTRUNGDANBENCHMARKS Test Trung's code on Dan's benchmarks
%   For single output, Trung's code is an improved version of Oppers
% This code (and Trung's code) uses the new version of gpml

DATADIR   = 'dataDan';
str_init  = datestr(now,30);
TRGDIR    = ['resultsTrung/', str_init];
mkdir(TRGDIR);
logFname  = [TRGDIR, '/', str_init, '.log']; 
diary(logFname);

%benchmarks = {'lineardata', 'poly3data', 'expdata', 'sindata', 'tanhdata'};
benchmarks = {'expdata', 'lineardata', 'poly3data',  'sindata', 'tanhdata'};

LEARN = 1; % loads results from file

for i = 1 : length(benchmarks)
  evaluateBenchmark(DATADIR, TRGDIR, benchmarks{i}, LEARN);
end
  

diary off;

return;

  
%% Evaluate a single benchmark
function evaluateBenchmark(DATADIR, TRGDIR, benchmark, learnFlag)
% Just avoids Matlab sending me stupid warning
f = [];  dfunc = []; noise = [];  train = [];
test = [];  x = [];  y = [];

        
load([DATADIR, '/', benchmark], 'f', 'dfunc', 'func', 'noise', 'train', ...
            'test', 'x', 'y');
train = train+1; test = test+1; % shifts indices to Matlab 
nFolds = size(train,1);
strFunc = parsePython2Matlab(func);
eval(['fwdFunc = @(theta)',strFunc, ';']);
for k = 1 :  nFolds
    [xtrain, ftrain, ytrain, xtest, ftest, ytest] = ...
        readSingleFold(x, f, y, train, test,k);
    fname = [TRGDIR, '/', benchmark, '_k', num2str(k), '.mat'];    
    if (learnFlag)
       [mufPred, sigmafStar, yStar, pred, m, conf] = runTrungSingle(fwdFunc, xtrain, ytrain, xtest, ...
                                                        noise);
        save(fname, 'mufPred', 'sigmafStar', 'yStar', 'pred', 'm', 'conf');
    else
        load(fname);
    end

    
  strTitle = [benchmark, '-Fold', num2str(k)];                                                  
  evaluateLatentPredictions(xtrain, ftrain, xtest, mufPred, sigmafStar, ...
                             ftest, strTitle);                                  
  evaluateObservablePredictions(xtrain, xtest, ytest, yStar, strTitle);   
  
end



%% Runs Trung's code on a Single Problem
%% Trung's code learns the noise so this parameter is not used
function [mufStar, sigmafStar, yStar, pred, m, conf] = runTrungSingle(fwdFunc, xtrain, ytrain, xtest,...
                                                        noise)

%% DELETE ME?
% ybar   = mean(ytrain);
% ytrain = ytrain - ybar; 
                                                    
% lets assume that we are given the noise for now
% We can have it as an additional parameter of the covariance
% May be not,  I think I need to learn it
[conf m] = getConfigTrung(fwdFunc, xtrain, ytrain, xtest);

%% Learn Posterior Parameters
 m                    = learnFullGaussian(m,conf);
[pred.mu, pred.sigma] = feval(m.pred, m, conf, m.xt); 
pred.Sigma            = diag(pred.sigma); % we don't really care about cross-covariances

%% DELETE ME?
% pred.mu = pred.mu + ybar;


mufStar               = pred.mu;
sigmafStar            = sqrt(diag(pred.Sigma));



%% We use my old function for predicting observables
conf.S          = conf.nsamples;
conf.elogMethod = 'mc';
yStar           = opperPredictObservable(pred, fwdFunc, {}, conf );
                                   
return;





function [conf m] = getConfigTrung(fwdFunc, xtrain, ytrain, xtest, noise)

covfunc = @covMatern52iso;

%% Configuration
conf.nsamples               = 10000;      % Number of samples for MC estimates
conf.covfunc                = covfunc;
conf.maxiter                = 200; % global iterations
conf.variter                = 500; % # iter var parameters
conf.hypiter                = 100; % # iter hyperparameters
conf.likiter                = 100; % # iter likelihood parameters
conf.displayInterval        = 20;
conf.checkVarianceReduction = false;
conf.learnhyp               = true;
conf.latentnoise            = 1e-10; % Variance of latent noise

%% Model 
myLogLike = @(y,f,hyp) llhGaussian(y,feval(fwdFunc,f),hyp);
[N D] =  size(xtrain); % D may be used by covfunc
m.x   = xtrain; m.y = ytrain; m.xt = xtest;
m.N   = N; 
m.Q   = 1; % Only a single latent function
                                          % initial  
m.pars.M = ytrain;                          % Variational

%m.pars.L = -2*(1./1e-3)*(ones(m.N*m.Q,1));  % parameters (-2*eta)
m.pars.L = -2*(1./1e10)*(ones(m.N*m.Q,1));

m.pars.hyp.covfunc = covfunc;   % covariance hyperparameters
m.pars.hyp.cov = cell(m.Q,1);
m.pars.hyp.cov{1}  = log(ones(eval(feval(covfunc)),1));
m.likfunc          = myLogLike;
m.pred             = @predRegression;


m.pars.hyp.likfunc = m.likfunc;
m.pars.hyp.lik     = 0.5*log(var(ytrain,1)/4 + 1e-4);
  


return;










  
