function [ nelbo, gradTheta ] = opperNelboGaussFull( theta,  prior, ...
                                                      fwdFunc, diagSigmay, y, S, method)
%NEGELBOOPPERGAUSSFULL Negative Evidence Lower Bound for Opper's method
%   Detailed explanation goes here
% S: The number of samples
% 
% theta = [nu; eta];
% mq     = (Sigmap)(omegaq) --> omegaq = (Sigmap{^-1})(mq)
% Sigmaq = (diag(eta) + So^{-1})^{-1} 
if (nargin == 7) 
    method = 'quad';
end

switch method,
        case 'quad',
            fnElog= @opperElogGaussFullQuad;
        case 'mc',
            fnElog= @opperElogGaussFullMC;
     otherwise,
        disp('Invalid ELOG Method');
end
[ omegaq, etaq ] = opperUnpack(theta);
post.omega = omegaq;
post.eta   = etaq;

%% Assign parameters from Prior
Sigmap   = prior.Sigma;
mup      = prior.mu;

%% Gets posterior parameters from the current representation
%[ muq, Sigmaq ] = opperGetMeanParam( omegaq, etaq, Sigmap, Lambdap );
[ muq, Sigmaq ] = opperGetMeanParamStruct( prior, post );
post.mu        = muq;
post.Sigma     = Sigmaq;


%% Evaluate the KL term
% kl  = klGauss( mup, Sigmap, muq, Sigmaq, 'chol' );
%kl   = klGaussStructChol(prior, post);
kl    = opperKLGaussReparam(prior, post);

%% Evaluate the likelihood term and its gradients if required
if (nargout == 1)
    elog = feval(fnElog, y, fwdFunc, muq, Sigmaq, diagSigmay, S);
else     
    [elog gradMuq, gradSigmaq] = feval(fnElog, y, fwdFunc, muq, Sigmaq, diagSigmay, S);
    
    % Computes the gradients of the Elbo
    gradOmegaq     = mup + Sigmap*(gradMuq - omegaq);
    gradEta     = (-0.5)*(Sigmaq.^2)*(etaq+2*gradSigmaq); 
    
    % Gradients of the negative elbo (nelbo)
    gradOmegaq    = - gradOmegaq;
    gradEta    = - gradEta;
    
    gradTheta  = opperGradientPack( gradOmegaq, gradEta, etaq );
end

%% Comptues elbo = variational lower bound
elbo = elog - kl;

nelbo = - elbo;




return;


