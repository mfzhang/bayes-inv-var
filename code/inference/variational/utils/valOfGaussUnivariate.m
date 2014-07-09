function  val = valOfGaussUnivariate( y, m, sigma2 )
%VALOFGAUSSUNIVARIATE Summary of this function goes here
%   Detailed explanation goes here
% y, m, sigma2 can be scalars or vectors


val =  (exp(-0.5*((y -m).^2)./sigma2))./sqrt((2*pi).*sigma2);



return;



