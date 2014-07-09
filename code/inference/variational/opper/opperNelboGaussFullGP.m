function [ nelbo gradLogHyper ] = opperNelboGaussFullGP( loghyper, meanFunc, ...
                                                    covFunc, fwdFunc, ...
                                                    prior, post, diagSigmay, ...
                                                    x, y, S, elogMethod)
%OPPERNELBOGAUSSFULLGP Nelbo for GP
%   Used Only for optimizing hyper-parameters
switch elogMethod,
        case 'quad',
            fnElog= @opperElogGaussFullQuad;
        case 'mc',
            fnElog= @opperElogGaussFullMC;
     otherwise,
        disp('Invalid ELOG Method');
end


elog = feval(fnElog, y, fwdFunc, post.mu, post.Sigma, diagSigmay, S);                          

   
if (nargin > 1)
[kl dKL] = opperKLGaussGPReparam( meanFunc, covFunc, loghyper, ...
                              x, post);                          
else
   kl  = opperKLGaussGPReparam( meanFunc, covFunc, loghyper, ...
                              x, post);                          
end

elbo  = elog - kl;
nelbo = -elbo;

gradLogHyper = dKL;



return;














