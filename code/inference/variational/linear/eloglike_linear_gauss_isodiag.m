function [ elog grad ] = eloglike_linear_gauss_isodiag(y, sigmay, ...
                                                    fwdfunc, fwd_param, ...
                                                    mu1, diagSigma1 )
%ELOGLIKE_LINEAR_GAUSS_ISODIAG: Expected loglikelihood for linearized model
% y: observations
% sigma_y: variance of likelihood model
% fwdfunc: Pointer to fwd function
% fwd_param: cell with parameters of fwd model
% mu1: posterior means
% diagSigma1: posterior variances
%
% requires DERIVESTsuite
%
N =size(y,1);

fobj = @(xx) fwdfunc(fwd_param, xx);

%% fwd model at posterior mean
[fwdval J] = feval(fobj, mu1); % second argument is the Jacobian

%% Estimates the Jacobian at the posterionr mean
%J    = jacobianest(fobj, mu1);

%% constant term
Z = (-N/2)*log(2*pi*sigmay);


%% quadratic term
lquad  = (y-fwdval)'*(y-fwdval);

%% trace term: trace(J'J*Sigmaq)
ltrace = sum(sum((J'*J).*diag(diagSigma1)));


elog = Z -  (1/(2*sigmay))*( lquad + ltrace);

end




