function [ y, param ] = generateDataisoLinear( N )
%GENERATEDATAISOLINEAR The output is a simple scaled version of the input
%   Detailed explanation goes here

if (nargin == 0)
    N = 100;
end

%% Prior parametes 
mup    = zeros(N,1);
Sigmap = rand(N);
Sigmap = Sigmap*Sigmap';

%% likelihood Parameters
fwdFunc  = @isoLinearFwdModel;
fwdParam =  10*rand(1,1);


%% Parameters of conditional likelihood
sigmay  = 1e-5;
lambday = 1/sigmay;
Sigmay  = sigmay*eye(N);


%% Draws one sample from M-dimensional prior
theta   = sampleGauss(mup,Sigmap,1); 


%% Draws from conditional likelihhood
f  = feval(fwdFunc, theta, fwdParam);
y  = sampleGauss(f, Sigmay, 1);
plot(theta, y, '.', 'MarkerSize', 12); 
xlabel('True Theta'); ylabel('Observations'); title('Data');


%% Get Natural Parameters
[ nup,  Lambdap, cholSigmap] = getMeanFromNaturalGauss( mup, Sigmap );


param.prior.nu        = nup;
param.prior.Lambda    = Lambdap;
param.prior.mu        = mup;
param.prior.Sigma     = Sigmap;
param.prior.cholSigma = cholSigmap;
%
param.like.lambda  = lambday;
param.fwd{1}       = fwdParam;
param.fwdFunc      = fwdFunc;

return;



