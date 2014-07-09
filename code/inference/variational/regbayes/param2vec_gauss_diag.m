function theta = param2vec_gauss_diag(param, grad_flag)
% reads cell of parameters into a vector
% the variances are given in log space
% for the variances: theta_i = log sigma_i --> theta_i = 0.5*log(Sigma_ii)
% grad_flag: flag indicating to retunrn the gradient of 
% the original parameter wrt the transformed parameter instead

if (nargin == 1 || grad_flag==0) % actual (transformed) parameters
    mu        = param{1};
    diagSigma = 0.5*log(param{2});
    theta = [mu(:); diagSigma(:)];

else % returns gradient of transformed parameter instead
    grad_mu        = ones(size(param{1}));
    grad_diagSigma = 2*param{2};
    theta = [grad_mu(:); grad_diagSigma(:)];
end
    


return;



