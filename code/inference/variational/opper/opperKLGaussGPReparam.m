function [ kl dKL gradK ] = opperKLGaussGPReparam( meanFunc, covFunc, loghyper, ...
                              x, post)
%OPPERKLGAUSSGPREPARAM KL Gauss GP using opper reparametrization
%   Detailed explanation goes here

N     = size(x,1);
mup   = meanFunc*ones(N,1); % Not really used (for now, mup = 0)
K     = feval(covFunc, loghyper, x);

%% Reads posterior parameters in Oppers parametrization
omegaq  = post.omega;
etaq    = post.eta;

%% Useful matrix for inverses
sqrtD   = diag(sqrt(etaq)); 
I       = eye(N);
E       = I + sqrtD*K*sqrtD;
cholE   = getChol(E);
Einv    = getInverseChol(cholE);


%% Trace Term
ltrace = trace(I - sqrtD*Einv*sqrtD*K);

%% Quadratic Term
lquad  = omegaq'*K*omegaq;

%% Determinant term
ldet   = getLogDetChol(cholE);

%% KL divergence
kl    = 0.5*(ltrace + lquad + ldet);

if (nargout > 1) % Gradients wrt hyper-parmeters are required
    nhyper = length(loghyper);
    gradK = cell(nhyper);
    Binv  = sqrtD*Einv*sqrtD;
    dKL = zeros(nhyper,1);
    for z = 1 : nhyper
        gradK{z}  = feval(covFunc, loghyper, x, z);
        dKL(z) =  0.5*trace((Binv - omegaq*omegaq')*gradK{z}); % TODO: More efficient
    end
end


return;


