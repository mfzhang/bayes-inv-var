function  smll = mySMLL(ytrain, ytest, muPred, varPred )
%MYSMLL Standardised mean log loss
%   ean Negative log likelihood aussuming gaussianity
%   and standardised by using a gaussian with the training data stats
mu     = mean(ytrain);
sigma2 = var(ytrain);

ll_model = - logOfGaussUnivariate( ytest, muPred, varPred );
ll_naive = - logOfGaussUnivariate( ytest, mu, sigma2 );

smll = mean(ll_model - ll_naive);


return;

