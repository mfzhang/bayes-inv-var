function [ nelbo ] = negElboLinearGaussGaussIso( y, fwdval, J, ...
                                                      mup, Sigmap, muq, Sigmaq, sigmay)
%NEGELBOLINEARGAUSSGAUSSISO Negative Evidence Lower Bound for linearized model
%   Detailed explanation goes here
% Negative Evidence Lower Bound for linearized model
% full Gaussian Posteriors and Isotropic likelihood

elog  = eLoglikeLinearGaussIsoFull(y, fwdval, J, muq, Sigmaq, sigmay);
kl  = klGauss( mup, Sigmap, muq, Sigmaq);

elbo = elog - kl;
nelbo = - elbo;


end

