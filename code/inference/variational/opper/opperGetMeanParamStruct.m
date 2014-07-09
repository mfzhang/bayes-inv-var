function [ mq, Sigmaq, Lambdaq ] = opperGetMeanParamStruct( prior, post )
%OPPERGETMEANPARAMSTRUCT It uses prior and posterior structues
%   


%% Get prior Parameters
K       = prior.Sigma;  % Prior Covariane
Lambdap = prior.Lambda; % Prior precision 
Deta    = diag(post.eta);
N       = length(post.eta);

%% Computes mean
mq            = K*post.omega; 


%% Posterior Precision
Lambdaq       = Lambdap + Deta;


%% 1st way to compute Sigmaq
% cholLambdaq   = getChol(Lambdaq);  
% Sigmaq        = getInverseChol(cholLambdaq);


%% 2nd way to computer Sigmaq
%Dinv = diag(1./post.eta);
%C      = (Dinv + K);
%cholC  = getChol(C);
%Cinv   = getInverseChol(cholC);
%Sigmaq        = K - K*Cinv*K; % (Kinv + Deta)^-1

%% 3rd way to computer Sigmaq
%Sigmaq  = Dinv*Cinv*K;

%% 4th way: hopefully most stable
sqrtD   = diag(sqrt(post.eta));
E       = eye(N) + sqrtD*K*sqrtD;
cholE   = getChol(E);
Einv    = getInverseChol(cholE);
Sigmaq  = K - K*sqrtD*Einv*sqrtD*K;



return;

