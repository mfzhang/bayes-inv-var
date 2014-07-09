function [ param, mu_q, sigma_q, nelbo ] = varLinearGaussStochastic( y, param, rbconf )
%VARLINEARGAUSSSTOCHASTIC Stochastic Optimization for Variational Linear
%Model

fprintf('using alpha=%f, nbatch=%d', rbconf.alpha, rbconf.nbatch);

%% Reads off  parameters of prior and likelihood
nu_p     = param.prior{1};
Lambda_p = param.prior{2};
lambday  = param.like{1};
[ mu_p, Sigma_p, cholL_p] = getMeanFromNaturalGauss(nu_p,Lambda_p);
D         = length(nu_p);
D2        = D.^2;
N         = length(y);

%% Initializes parameters and assings best so far
[nu_q0 Lambda_q0 cholL_q0 mu_q0 ] = initVariationalGaussNatural(D,param.post);
nu_q      = nu_q0;
mu_q      = mu_q0;
Lambda_q  = Lambda_q0;
cholL_q   = cholL_q0;


%% Here we iterate over mini-batches
% The idea is that these operations can be distributed 
i   = 0;
err = inf;
while ((i < rbconf.maxiter) && (err > rbconf.tol))
    mu_qOld      = mu_q;
    Lambda_qOld  = Lambda_q;
    %nu_qOld      = nu_q;
    
    % Selects subset of observations (can do this outside the loop)
    idx      = randperm(N, rbconf.nbatch)';
    ybatch   = y(idx);
    
    %% Evaluates fwd model and returns evaluations at batch
    [ fwdval, J] = wrapNonlinearFwdBatch( rbconf.fwdfunc, mu_q, idx );
 
    %% Noisy gradient wrt precision and updates
    gradLambda = Lambda_p - Lambda_q + lambday*(J'*J);
    Lambda_q   = Lambda_q + rbconf.alpha*gradLambda;
    
    %% Noisy gradient wrt precision-adjusted mean and updates
    gradNu =  nu_p - Lambda_p*mu_q + lambday*J'*(ybatch - fwdval);
    nu_q   = nu_q + rbconf.alpha*gradNu;
    
    %% Updates mu_q
    cholL_q   = chol(Lambda_q, 'lower');
    mu_q      = solve_chol(cholL_q',nu_q);

    errMu     =  sum((mu_q - mu_qOld).^2)/D;    
    errLambda = sum( sum((Lambda_q - Lambda_qOld).^2) )/D2;
    err       =  max([errLambda, errMu]);
    %err = errMu;
    
    i = i+1;
    showProgressVariational(i, [], err);
end


%% For analysis purposes we compute SigmaPred
Sigma_q = cholL_q'\(cholL_q\eye(D));
sigma_q = diag(Sigma_q);

%% Assign parameters to output
param.post = {nu_q, Lambda_q};

nelbo = [];

return;




