function [ kl gradMuq gradSigmaq ] = evbKlGaussFull( mup, Sigmap, muq, Sigmaq)
%EVBKLGAUSSFULL Computes KL divergence KL(q | p) and its gradients
%   for full Gaussian case

%% KL divergence
[kl Lambdap, Lq]  = klGauss( mup, Sigmap, muq, Sigmaq, 'chol' );
        

%% gradient wrt Muq
gradMuq = -Lambdap*(mup - muq);

%% Gradient wrt Sigmaq
Lambdaq  = getInverseChol(Lq);
gradSigmaq = (0.5)*(Lambdap - Lambdaq);  



return;

