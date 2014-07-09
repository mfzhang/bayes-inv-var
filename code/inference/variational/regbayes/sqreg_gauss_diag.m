function [lreg grad] = sqreg_gauss_diag(ptr_fwd, fwd_param, ...
                                 A, sigmay, mu1, diagSigma1, lambda_reg)
% Squared regularizer for approximate Gaussian likelihood model
% p(y | ..) = N(y | A*theta, sigmay)
% the forward model is evaluated in the space of parameters theta
% here fwd_param indicates additional (fixed) parameters to the forward
% model (e.g. sensitivity matrix, etc)
% lambda_reg is the regularization parameter
% mu1:       Current mean of posterior
% diagSigma1: current (diagonal) posteriot covaraindce 
% 28/10/2013
% lambda_reg is used outside this function

fwdval = feval(ptr_fwd, fwd_param, mu1);
lreg   = sum((A*mu1 - fwdval).^2);

if (nargout == 2) % returns derivatives here wrt posterior parameters
    grad = zeros(numel(mu1)+numel(diagSigma1),1);
end

%% DELETE ME
%lreg = 0;

return;

