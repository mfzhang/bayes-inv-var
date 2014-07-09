function [elog grad] = eloglike_gauss_iso_diag(y, A, sigmay, mu1, diagSigma1)
% expected log likelihood when conditional likelihood is isotropic
% and posterior is diagonal gaussian
% the likelihood is "approximated with a linear model A*theta
%
% y: is a column vector
% sigmay: variance of observations

N = size(y,1);

%% Precompute for convenience
%C = A'*A;
alpha_mu = A*mu1;
prod_mu  = A'*alpha_mu;
diagC    = diagProd(A',A);

%% constant term
Z = (-N/2)*log(2*pi*sigmay);

%% quadratic term
lquad = y'*y - 2*y'*alpha_mu + mu1'*prod_mu;

%% trace term
ltr = diagC'*diagSigma1;

elog = Z - (1/(2*sigmay))*(lquad + ltr);


%% if required it returns the gradients
if (nargout  == 2) % returns gradient if required
    dl_dm      = (1/sigmay)*(A'*y - prod_mu);
    dl_dS      = (-0.5/sigmay)*diagC;
    grad_l     = [dl_dm(:); dl_dS(:)];
    grad_theta = param2vec_gauss_diag({mu1, diagSigma1}, 1); % ask for gradient of parametes
    grad       = grad_l.*grad_theta;
end

return;


