function [kl grad] = kl_gauss_iso_diag(mu0, sigma0, mu1, diagSigma1)
% KL divergence when prior is isotropic Gaussian and posterior 
% is diagonal Gaussian
% kl(p(mu0, sigma^2)| q(mu1, diag(diagSima1) )
% mu0: mean of p
% Sigma0: variance of p
% mu1: mean of q
% diagSigma1: column vector of q variances

D = length(mu0);


%% Right KL
%
invsigma0 = 1/sigma0;

%% trace term
ltr = invsigma0*sum(diagSigma1);

%% Quadratic term
lquad = sum(((mu1 - mu0).^2)/sigma0);

%% Log determinant term
ldet = D*log(1/sigma0) + sum(log(diagSigma1));

%% final result
kl = 0.5*(ltr + lquad - ldet - D);

%% returns gradients if required
if (nargout == 2)
    dkl_dm = invsigma0*(mu1 - mu0);
    
    %% 
    invdiagSigma1 = 1./diagSigma1;
    dkl_dS  = 0.5*( invsigma0 - invdiagSigma1 );
    
    grad_l     = [dkl_dm(:); dkl_dS(:)];
    grad_theta = param2vec_gauss_diag({mu1, diagSigma1}, 1); 
    grad       = grad_l.*grad_theta;
end



%% Below the wrong KL!
% invdiagSigma1 = 1./diagSigma1; % vector with the inverse (diagonal)
% 
% %% trace term
% ltr = sigma0*sum(invdiagSigma1);
% 
% %% quadratic term
% lquad = sum(((mu1 - mu0).^2)./diagSigma1);
% 
% %% log determinant term
% ldet = D*log(sigma0) + sum(log(invdiagSigma1));
% 
% %% final result
% kl = 0.5*(ltr + lquad - ldet - D);
% 
% %% returns gradients if required
% if (nargout == 2)
%     dkl_dm = invdiagSigma1.*(mu1 - mu0);
%     
%     invd2   = invdiagSigma1.^2; % squared of inverse
%     dkl_dS  = 0.5*( -sigma0*invd2 -  ((mu1 - mu0).^2).*invd2 ...
%                 + invdiagSigma1 );
%     
%     grad_l     = [dkl_dm(:); dkl_dS(:)];
%     grad_theta = param2vec_gauss_diag({mu1, diagSigma1}, 1); 
%     grad       = grad_l.*grad_theta;
% end

return;
