function  Perf = evaluateOpperDan(  )
%EVALUATEOPPERDAN  Evaluate predictions on Dan's benchmarks
%   Detailed explanation goes here
DATADIR    = 'dataDan';
%RESULTSDIR =  'resultsTrung/20141004T153125';
% RESULTSDIR =  'resultsTrung/20141002T172455';
RESULTSDIR =  'resultsTrung/20141010T182138';

%RESULTSDIR =  'resultsTrung/old';
%RESULTSDIR =  'resultsOpper';


benchmarks = {'lineardata', 'poly3data', 'expdata', 'sindata', 'tanhdata'};
%benchmarks = {'lineardata', 'poly3data', 'expdata'};
%benchmarks = {'lineardata'};



L = length(benchmarks);

Perf = zeros(L,6);
fprintf('Benchmark\tNLL-MEAN\tNLL-STD\tSMSE-MEAN-F\tSMSE-STD-F\tSMSE-MEAN-Y\tSMSE-STD-Y');
fprintf('\n');
for i = 1 : L
  [Perf(i,1), Perf(i,2), Perf(i,3), Perf(i,4) Perf(i,5) Perf(i,6)] = evaluateBenchmark(DATADIR, benchmarks{i}, RESULTSDIR);

  fprintf('[9] %s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f',benchmarks{i}, Perf(i,1), Perf(i,2), Perf(i,3), Perf(i,4), Perf(i,5), Perf(i,6) );
  fprintf('\n');
end
 

return;






%% Single benchmark
function [meanMLL stdMLL meanSMSE stdSMSE  meanSMSEy stdSMSEy] = evaluateBenchmark(DATADIR, benchmark, RESULTSDIR)
% Just avoids Matlab sending me stupid warning
f = [];  dfunc = []; noise = [];  train = [];
test = [];  x = [];  y = [];

        
load([DATADIR, '/', benchmark], 'f', 'dfunc', 'func', 'noise', 'train', ...
            'test', 'x', 'y');
train = train+1; test = test+1; % shifts indices to Matlab 
[nFolds, nTrain] = size(train);
vSMSE  = zeros(nFolds,1);
vSMSEy = zeros(nFolds,1); % ystar
vMLL   = zeros(nFolds,1);
for k = 1 : nFolds
    % Get true test data
    [xtrain, ftrain, ytrain, xtest, ftest, ytest] = ...
               readSingleFold(x, f, y, train, test,k);
    
  % Get predictions
  fname = [RESULTSDIR, '/', benchmark, '_k', num2str(k), '.mat'];
  load(fname, 'mufPred', 'sigmafStar', 'yStar', 'pred');
  
  % Evaluate performance
  vSMSE(k)  = mySMSE(ftrain, ftest, mufPred);
  vMLL(k)   = myMLL(ftrain, ftest, pred.mu, diag(pred.Sigma));
  vSMSEy(k) = mySMSE(ytrain, ytest, yStar);  
end

meanSMSE = mean(vSMSE);
stdSMSE  = std(vSMSE);
meanMLL  = mean(vMLL);
stdMLL   = std(vMLL);
meanSMSEy = mean(vSMSEy);
stdSMSEy  = std(vSMSEy);




  
return;
  
