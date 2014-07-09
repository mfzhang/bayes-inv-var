function [ elog ] = eLoglikeLinearGaussIsoFull(y, fwdval, J, muq, Sigmaq, sigmay)
%ELOGLIKELINEARGAUSSISOFULL Summary of this function goes here
%   Detailed explanation goes here
% Explected log likelihood of a Linearized posterior 
% for an Isotropic Gaussian Likelihood centered at the linearized
% forward model and a Full posterior
% y: Actual observations
% fwdval: fwd model evaluated at muq
% J: Jacobian of fwd model evalauted at muq
% muq:  mean and covariance of posterior (not used here)
% Sigmaq:  covariance of posterior
% sigmay: noise variance
N =size(y,1);

%% constant term
Z = - (N/2)*log(2*pi*sigmay);

%% quadratic term
lquad  =  - (1/(2*sigmay))*(y-fwdval)'*(y-fwdval);

%% trace term: trace(J'J*Sigmaq)
ltrace = -  (1/(2*sigmay))*sum(sum((J'*J).*Sigmaq));

elog = Z +  lquad  + ltrace;


end

