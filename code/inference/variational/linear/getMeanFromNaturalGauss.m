function [ mu, Sigma, cholLambda] = getMeanFromNaturalGauss( nu, Lambda )
%GETMEANFROMNATURAL Gets the mean Parameters from Natural Parameters
%   nu: Precision-adjusted mean
%   Lambda: Precision
%
%   mu: Mean
%   Sigma: Covariance
%   cholLambda: Cholesky decomposition of Lambda
% We can use it to get the natural from the mean parameters as well

if (isDiag(Lambda))
    cholLambda     = sqrt(Lambda);
    Sigma = diag(1./diag(Lambda));
    mu    = Sigma*nu;
else
    cholLambda      = chol(Lambda, 'lower');
    mu              = solve_chol(cholLambda',nu);
    Sigma           = getInverseChol(cholLambda);
end

return;






