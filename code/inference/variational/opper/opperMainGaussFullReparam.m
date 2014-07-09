function [ param, muq, sigmaPred, nelbo ] = opperMainGaussFullReparam( y, param, conf )
%OPPERMAINGAUSSFULLREPARAM Summary of this function goes here
%   This one uses the reparametrization directly
%
optimizer = 'minFunc';

%% Read first configuration
S          = conf.S;
fwdFunc    = @(xx) param.fwdFunc(xx, param.fwd{:});
elogMethod = conf.elogMethod;


%% Reads  likelihood parameters
N       = size(y,1);
lambday    = param.like.lambda;
sigmay     = 1/lambday;
diagSigmay = sigmay*ones(N,1);


 %% Parameter initialization
[omega0 eta0]         = opperInitParametersReparam(param.post, N);
theta0                = opperPack( omega0, eta0 );

%% Optimization objective
fobj = @(theta) opperNelboGaussFull( theta,  param.prior, ...
                                    fwdFunc, diagSigmay, y, S, elogMethod);
                                




%% Optimization here
optconf = opperGetOptimizerConf(optimizer, conf);
[theta, nelbo] = optimizeObjective(optimizer, fobj, theta0, optconf);


%% Gets mean parameters
[ omegaq, etaq ]  = opperUnpack( theta ); %% this omegaq is Opper parameterization 
param.post.omega  = omegaq;
param.post.eta    = etaq;



%[ muq, Sigmaq Lambdaq ] = opperGetMeanParam( omegaq, etaq, Sigmap, Lambdap );
[ muq, Sigmaq Lambdaq ] = opperGetMeanParamStruct( param.prior, param.post);

sigmaPred         = diag(Sigmaq);
param.post.nu     = Lambdaq*muq; 
param.post.Lambda = Lambdaq;
param.post.mu     = muq;
param.post.Sigma  = Sigmaq;


return;


  