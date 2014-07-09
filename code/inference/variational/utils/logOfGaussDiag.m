function logn = logOfGaussDiag( y, m, diagS )
%LOGOFNORMALISO Log of Normal for Diagonal covariance case
% y: observations
% m: mean
% diagS: diagonal vector with variances
N = size(y,1);

logn = -0.5*(N*log(2*pi) + sum(log(diagS)) + sum(((y-m).^2)./diagS) );




return;

