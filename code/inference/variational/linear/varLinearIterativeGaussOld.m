function [ param, mu_q, sigmaPred, nelbo err] = varLinearIterativeGaussOld( y, param, rbconf)
%VARLINEARITERATIVEGAUSS Summary of this function goes here
% all the parametrization for this function are in CANONICAL form
% i.e. nu and Precision
%   Closed-form iterative updates for Gaussian Posterior

SHOW_NELBO = 0;


%% Reads off  parameters of prior and likelihood
nu_p     = param.prior{1};
Lambda_p = param.prior{2};
lambday  = param.like{1};
%[ mu_p, Sigma_p, cholL_p] = getMeanFromNaturalGauss(nu_p,Lambda_p);

%% Pointer to fwd function
fobj  = @(xx) rbconf.fwdfunc(xx, param.fwd{:});

%% Initializes parameters 
[nu_q0 Lambda_q0 cholL_q0 mu_q0 ] = initVariationalGaussNatural(rbconf.D,param.post);
nu_q      = nu_q0;
mu_q      = mu_q0;
Lambda_q  = Lambda_q0;
cholL_q   = cholL_q0;


%% Main Loop from now
i     = 1;
nelbo = NaN*ones(rbconf.maxiter,1);
err   = NaN*ones(rbconf.maxiter,1);
err(1)= inf;
%
% Variational objective at initialization
%nelbo(i) = getVariationalObjective(y, fobj, nu_p, Lambda_p, mu_q, cholL_q, lambday);
showProgressVariational(i, nelbo(i), err(i));
%min_nelbo = nelbo(i);
%
while ((i < rbconf.maxiter) && (err(i) > rbconf.tol))
    muOld      = mu_q;
    
    %% Fwd models and its Jacobian wrt mu
    [fwdval J] = feval(fobj, mu_q);
    
    %% Update nu
    nu_q     = nu_p + J'*lambday*(y - fwdval + J*mu_q);
        
    %% Update Precision
    Lambda_q = lambday*(J'*J) + Lambda_p;

    
    %% Recompute mean parameters (mu need to evaluate fwd model)
    cholL_q   = chol(Lambda_q, 'lower');
    mu_q      = solve_chol(cholL_q',nu_q);
          
    %% Looks at Variational Objective
    i = i+1;
    % nelbo(i) = getVariationalObjective(y, fobj, nu_p, Lambda_p, mu_q, cholL_q, lambday);
    err(i)   = sum((mu_q - muOld).^2)/rbconf.D;
    showProgressVariational(i, nelbo(i), err(i));
    
%     if ( nelbo(i) < min_nelbo) % improved variational objective
%         fprintf('*** Variational Objective improved ***\n');
%         nuBest      = nu_q;        
%         muBest      = mu_q;
%         LambdaBest  = Lambda_q;
%         min_nelbo    = nelbo(i);
%     end
    
%     if ( mod(i,10) == 0 )
%         save(['tmp_',num2str(i), '_',datestr(now),'.mat']);
%     end
        
end

%% Update final output with  best values
% nu_q      =  nuBest;        
% mu_q      =  muBest;
% Lambda_q  =  LambdaBest;

        
%% For analysis purposes we compute SigmaPred
%cholL_q  = chol(Lambda_q, 'lower');
SigmaPred = cholL_q'\(cholL_q\eye(rbconf.D));
sigmaPred = diag(SigmaPred);

%% Assign parameters to output
param.post = {nu_q, Lambda_q};


%% Gets rid off "extra" iterations
%nelbo(isnan(nelbo)) = [];
err(isnan(err)) = [];


return;





