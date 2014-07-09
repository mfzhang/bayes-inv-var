function [ kl dKL] = klGaussGP(meanFunc, covFunc, loghyper, x,  muq, Sigmaq )
%KLGAUSSGP KL(q|p) where q is a Guassian and p is a GP prior
N = size(x,1);
strDecomp = 'chol';

mup = meanFunc*ones(N,1);
K   = feval(covFunc, loghyper, x);
if (nargout > 1) % Gests the fradients of the hyperparameters
    [kl invK] = klGauss( mup, K, muq, Sigmaq, strDecomp );
    nhyper = length(loghyper);
    dKL  = zeros(nhyper,1);
    for z = 1 : nhyper
        gradK = feval(covFunc, loghyper, x, z);
        v     = invK*(mup-muq); 
        dKL(z)   = 0.5*( ...
                        -traceProd(invK*Sigmaq,invK*gradK) ...
                        - v'*gradK*v ...
                        + traceProd(invK, gradK) ...
                        );
    end
else % only KL is needed
    kl = klGauss( mup, K, muq, Sigmaq, strDecomp );
end
    
    
return;





