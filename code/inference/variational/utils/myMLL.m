function  mll  = myMLL( ytrain, ytest, muPred, varPred )
%MYMLL Mean Log Loss 
%   Mean Negative log likelihood aussuming gaussianity


ll = - logOfGaussUnivariate( ytest, muPred, varPred );

mll = mean(ll,1);

return;






