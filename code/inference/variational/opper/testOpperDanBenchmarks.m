function testOpperDanBenchmarks()
% Test oppers alg. with Dan's bencgmarks
clear all; clc; close all;
DATADIR = 'dataDan';
benchmarks = {'lineardata'};


for i = 1 : length(benchmarks)
  evaluateBenchmark(DATADIR, benchmarks{i});

end
 
return;



%% 
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
    
   [mufPred, sigmafStar, yStar, pred] = runOpperSingle(fwdFunc, xtrain, ytrain, xtest, ...
                                                        noise);
    
   evaluateLatentPredictions(xtrain, ftrain, xtest, mufPred, sigmafStar, ...
                             ftest);                                  
                         
  evaluateObservablePredictions(xtrain, xtest, ytest, yStar);   
  
  fname = ['resultsOpper/', benchmark, '_k', num2str(k), '.mat'];
  save(fname, 'mufPred', 'sigmafStar', 'yStar', 'pred');
end


  
return;
  


%% Read data
function [xtrain, ftrain, ytrain, xtest, ftest, ytest] = ...
               readSingleFold(x, f, y, train, test, k)
    xtrain = x(train(k,:),:);
    ftrain = f(train(k,:),:);
    ytrain = y(train(k,:),:);

    xtest = x(test(k,:),:);
    ftest = f(test(k,:),:);
    ytest = y(test(k,:),:);
    
return;


%% Evaluate Latent Predicitons
function evaluateLatentPredictions(xtrain, ftrain, xtest, muStar, sigmaStar, ftest)
figure;
if size(xtest,1)
    [foo idx] = sort(xtest);
    plot_confidence_interval(xtest(idx),muStar(idx),sigmaStar(idx),1.96);
    hold on;
    plot(xtest(idx), ftest(idx), 'g--', 'MarkerSize', 12);
    plot(xtrain, ftrain, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    xlabel('x');
    ylabel('f');
    title('Latent Function');
    legend({'Pred Error Bars', 'Pred Mean', 'True Test', 'Training'});
    
else
    disp('This problem is not 1-dimesional, plots not shown');
end

return;


%% 
function evaluateObservablePredictions(xtrain, xtest, ytest, yStar)
if (size(xtrain,2) == 1)
    [foo idx] = sort(xtest);
    xtest = xtest(idx);
    ytest = ytest(idx);
    yStar = yStar(idx);
    figure;
    plot(xtest, yStar, 'k-', 'LineWidth', 2, 'MarkerSize', 12);
    hold on;
    plot(xtest, ytest, 'g--', 'MarkerSize', 12);
    hold on;
    plot(xtrain, min(ytest)*ones(size(xtrain,1),1), 'rx', 'MarkerSize', 12);

    title('Observable Function');
    legend({'Predicted', 'True', 'Training Input Location'});
end

return;







%% Runs Opper on a Single Problem
function [mufStar, sigmafStar, yStar, pred] = runOpperSingle(fwdFunc, xtrain, ytrain, xtest,...
                                                        noise)
% lets assume that we are given the noise for now
% We can have it as an additional parameter of the covariance
% May be not,  I think I need to learn it
[conf param] = getConfigOpper(fwdFunc, noise);

%% Learn Posterior Parameters
 param  = opperMainAll(xtrain,ytrain, param, conf);
pred   = opperPredictLatent( param.meanFunc, param.covFunc, param.loghyper, ...
                           param.prior, param.post, xtrain, xtest);                                     
mufStar    = pred.mu;
sigmafStar = sqrt(diag(pred.Sigma));

yStar     = opperPredictObservable(pred, param.fwdFunc, param.fwd, conf );
                                   
return;




%%
function [conf param]  = getConfigOpper(fwdFunc, std_noise)
conf.MaxIter       = 1000;
conf.alpha         = 0.01; % learning rate
conf.S             = 10000; % Number of samples for MC estimate
conf.elogMethod    = 'mc'; % 'quad' for quadrature or 'mc' for montecarlo

param.meanFunc    = 0;
% param.covFunc    = 'covSEard';
param.covFunc     = 'covMatern5iso';
param.fwd         = {};
param.fwdFunc     = fwdFunc;
lambday           = 1/(std_noise^2);
param.like.lambda = lambday;
param.post        = {};


return;




