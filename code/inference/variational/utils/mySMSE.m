function  e  = mySMSE(ytrain, ytest, ypred )
% function r = mse(ytest,ypred)
% computes the SMSE: The mean square error normalised by the variance
% of the test targets
% INPUT:
%   - ytest: real
%   - ypred: prediction from my model
%
% Edwin V. Bonilla

res = ytest - ypred; % residuals
e   = mean(res.^2,1);

flagNormalization = 1;                   % 0: uses N-1;   1: uses N
vari = var(ytest, flagNormalization, 1); % Normalizes over N
e = e./vari;

return;



