function logval  = logOfGaussUnivariate( y, m, sigma2 )
%LOGOFGAUSSUNIVARIATE Summary of this function goes here
%   Detailed explanation goes here

logval = -0.5*(log(2*pi*sigma2) + ((y -m).^2)./sigma2);

return;

