function  [omega, eta]  = opperInitParameters( Lambdap, paramPost )
%OPPERINITETA Initializes vector that fully parametrize posterior
%covariance
% D: Dimensionality of covarince
% Lambda_q = (Lambda_p + diag(eta))
% mq   = Sigmap*nuq --> nuq = Lambdap*mq 
% The parametrization useD:

D = size(Lambdap,1);

if (isempty(paramPost))
    eta    = zeros(D,1);           
    mu     = randn(D,1); % Posterior mean
    omega  = Lambdap*mu;
else % simply assigns whatever is given
    omega  = Lambdap*paramPost.mu; 
    eta    = diag(paramPost.Lambda - Lambdap);
end



% I tried this before but ...
% eta                   = diag(Lambdaq - Lambdap);


% another idea:
% reasoning: eta is a small perturbarcion on the precisions
%var = 1e-7;
%eta = (1/var)*ones(D,1);



return;
