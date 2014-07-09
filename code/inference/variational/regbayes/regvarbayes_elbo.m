function [elbo grad] = regvarbayes_elbo(y, ptr_kl, ptr_elike, ptr_reg, ...
                                        ptr_fwd, prior_param, like_param, ...
                                        post_param, fwd_param, reg_param )
% Regularized Variational Bayes evidence lower bound
% ptr_kl: pointer to KL(q||p)
% ptr_like: pointer to expected log likelihood
% ptr_reg: pointer to regularizer
% ptr_fwd: pointer to fwd model
% prior_param: cell of prior parameters
% post_param: cell of posterior parameters
% y: observations
% fwd_param: parameters of fwd model
% It assumes we are approximating the fwd model with a lineat A*theta
% where theta are the parameters of interest
% and the actual forward model is used as a regularizer

lambda_reg = reg_param{1};


kl_param    = [prior_param post_param];
elog_param  = [like_param, post_param]; 
reg_param   = [fwd_param, like_param, post_param, reg_param];
    
if (nargout == 2) % gets the gradient
    [kl grad_kl]     = feval(ptr_kl, kl_param{:});
    [elog grad_elog] = feval(ptr_elike, y, elog_param{:});
    [reg grad_reg]   = feval(ptr_reg, ptr_fwd, reg_param{:});

    % modified on 28/10/2013
    grad             = grad_elog - grad_kl - (lambda_reg)*grad_reg;
    % Convex formulation
    %grad             = (1-lambda_reg)*(grad_elog - grad_kl) - (lambda_reg)*grad_reg;
else    
    kl          = feval(ptr_kl, kl_param{:});    
    elog        = feval(ptr_elike, y, elog_param{:});    
    reg         = feval(ptr_reg, ptr_fwd, reg_param{:});
end

% modified on 28/10/2013
elbo  = elog - kl - (lambda_reg)*reg;
% Convex formulation
%elbo       = (1-lambda_reg)*( elog - kl) - (lambda_reg)*reg;

return;

