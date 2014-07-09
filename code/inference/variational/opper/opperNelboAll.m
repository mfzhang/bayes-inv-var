function [ nelbo gradTheta ] = opperNelboAll( theta, meanFunc, covFunc, fwdFunc, ...
                                          diagSigmay, x, y, S, elogMethod)
%OPPERNELBOALL Computes the Nelbo using all parameters
%theta      = [omega; log(eta); loghyper];
N          = size(x,1);    % Number of data-points

idxOmega   = 1:N;
idxEta     = N+1:2*N;
idxVar     = [idxOmega, idxEta];
idxHyper   = 2*N+1:length(theta);

post_param = theta(idxVar);
loghyper   = theta(idxHyper);


prior.mu = meanFunc*ones(N,1); % prior mean
prior    = opperRecomputePrior(covFunc, loghyper, prior, x);

%% Fills in posterior
[ post.omega, post.eta ] = opperUnpack( post_param );
[ post.mu, post.Sigma ]   = opperGetMeanParamStruct( prior, post );


switch elogMethod,
        case 'quad',
            fnElog= @opperElogGaussFullQuad;
        case 'mc',
            fnElog= @opperElogGaussFullMC;
     otherwise,
        disp('Invalid ELOG Method');
end



if (nargin > 1)
   [elog gradMuq, gradSigmaq] = feval(fnElog, y, fwdFunc, post.mu, post.Sigma, diagSigmay, S);
   [kl, dKL, gradK]   = opperKLGaussGPReparam( meanFunc, covFunc, loghyper, ...
                              x, post);     
                          
    % Computes the gradients of the Elbo wrt variational parameters
    gradOmegaq     = prior.mu + prior.Sigma*(gradMuq - post.omega);
    gradEta     = (-0.5)*(post.Sigma.^2)*(post.eta + 2*gradSigmaq); 
    
    % Gradients of the negative elbo (nelbo)
    gradOmegaq    = - gradOmegaq;
    gradEta       = - gradEta;    
    gradVar       = opperGradientPack( gradOmegaq, gradEta, post.eta );    
    
    %% assing corrected gradients to be used below for implicit grads
    % Actually we shouldn't do this 
    %gradOmegaq       = gradVar(idxOmega);
    %gradEta          = gradVar(idxEta);
   
    gradHyper     = dKL;   
%     %% Additional implict gradients of hyper from nelog
%     for z = 1 : length(gradK)
%         gradHyper(z) = gradHyper(z)  ...
%                         + gradOmegaq'*(prior.Lambda*gradK{z}*post.omega) ...
%                         + gradEta'*diagProd(prior.Lambda*gradK{z},prior.Lambda);
%     end
     gradTheta = [gradVar;gradHyper];
   
else
    elog = feval(fnElog, y, fwdFunc, post.mu, post.Sigma, diagSigmay, S);   
    kl   = opperKLGaussGPReparam( meanFunc, covFunc, loghyper, ...
                              x, post);  
end


elbo  = elog - kl;
nelbo = - elbo;



return;



