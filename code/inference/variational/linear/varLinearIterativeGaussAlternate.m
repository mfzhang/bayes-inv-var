function [ param, mu_q, sigma_q, nelbo ] = varLinearIterativeGaussAlternate( y, param, rbconf)
%VARLINEARITERATIVEGAUSSALTERNATE 
% all the parametrization for this function are in CANONICAL form
% i.e. nu and Precision
%   Closed-form iterative updates for Gaussian Posterior
% Alternate optimization of parameters

%% Reads off  parameters of prior and likelihood
nu_p     = param.prior{1};
Lambda_p = param.prior{2};
lambday  = param.like{1};
[ mu_p, Sigma_p, cholL_p] = getMeanFromNaturalGauss(nu_p,Lambda_p);
D         = length(nu_p); 

%% Pointer to fwd function
%fobj  = @(xx) rbconf.fwdfunc(xx);
fobj = @(xx) wrapNonlinearFwd(rbconf.fwdfunc, xx);
%rbconf.fwdfunc = @(xx) wrapNonlinearFwdAdiff(theta2Grav, xx);

%% Initializes parameters and assings best so far
[nu_q0 Lambda_q0 cholL_q0 mu_q0 ] = initVariationalGaussNatural(D,param.post);
nu_q      = nu_q0;
mu_q      = mu_q0;
Lambda_q  = Lambda_q0;
cholL_q   = cholL_q0;


%% Main Loop from now
i     = 1;
nelbo = NaN*ones(rbconf.maxiter,1);
err   = inf;
%
% Variational objective at initialization
[fwdval J] = feval(fobj, mu_q); % fwd model and Jacobian
nelbo(i)   = getVariationalObjectiveFixed(y, fwdval, J, ...
                nu_p, Lambda_p, mu_q, cholL_q, lambday);
showProgressVariational(i, nelbo(i), err);
while ((i < rbconf.maxiter) && (err > rbconf.tol))
    
    %% Update Precision and its Cholesky
    Lambda_q = lambday*(J'*J) + Lambda_p;
    cholL_q   = chol(Lambda_q, 'lower');
    nl    = getVariationalObjectiveFixed(y, fwdval, J, ...
                                 nu_p, Lambda_p, mu_q, cholL_q, lambday );
    fprintf('Precision updated \n');
    showProgressVariational(i, nl, abs(nl-nelbo(i)));
   
    %% Update Precisio-adjusted mean
    % --> mu_q = mu_p + Sigma_p*J'*lambday*(y - fwdval); <-- DOESN"T WORK!
    %nu_q       = nu_p + J'*lambday*(y - fwdval + J*mu_q);
    %mu_q       = solve_chol(cholL_q',nu_q);
    mu_q = optimizeMeanVarGauss(y, fobj, mu_p, cholL_p, mu_q, cholL_q, lambday);
    fprintf('Mean updated \n');                             

   %% Update fwd model and Jacobian at the new mean 
   [fwdval J] = feval(fobj, mu_q); % fwd model and Jacobian

    
    %% Show progress
    i = i + 1;
    nelbo(i) = getVariationalObjectiveFixed(y, fwdval, J, ...
                                 nu_p, Lambda_p, mu_q, cholL_q, lambday );    
    err   = abs(nelbo(i) - nelbo(i-1));
    
    showProgressVariational(i, nelbo(i), err);
    
    if ( mod(i,10) == 0 )
        save(['tmp/tmp_',num2str(i), '_',datestr(now),'.mat']);
    end
    
end


%% Computes precision
nu_q = Lambda_q*mu_q;


%% For analysis purposes we compute SigmaPred
Sigma_q = cholL_q'\(cholL_q\eye(D));
sigma_q = diag(Sigma_q);

%% Assign parameters to output
param.post = {nu_q, Lambda_q};


%% Gets rid off "extra" iteration
nelbo(isnan(nelbo)) = [];

return;







