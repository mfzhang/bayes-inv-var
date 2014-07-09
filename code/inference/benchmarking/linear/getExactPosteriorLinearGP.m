function [ muq, Sigmaq ] = getExactPosteriorLinearGP( param, x , y )
%GETEXACTPOSTERIORLINEARGP Summary of this function goes here
%   Detailed explanation goes here

N = size(y,1);

%% prior
nup       = param.meanFunc*ones(N,1);
Sigmap    = feval(param.covFunc, param.loghyper, x);
cholSigmap = getCholSafe(Sigmap);
Lambdap   = getInverseChol(cholSigmap);

[mup, Sigmap] = getMeanFromNaturalGauss(nup, Lambdap);

%% Likelihood
A       = param.fwd{1};
[cols rows]   = size(A);
if (cols> 1 && rows==1) % a vector has been specified
    A = diag(A);
end

lambday = param.like.lambda;
sigmay  = 1/lambday;
Sigmay  = sigmay*eye(N);
lambday = lambday*eye(N); 

%% Posterior mean
%C     = A*Sigmap*A' + Sigmay;
%cholC = getChol(C);
%mu    = mup + (A*Sigmap)'*solveLinChol(cholC,y); 
Lambdaq     = Lambdap + A'*lambday*A;
cholLambdaq = getChol(Lambdaq);
Sigmaq      = getInverseChol(cholLambdaq);
muq         = Sigmaq*(A'*lambday*y + Lambdap*mup);

%% Posterior covariance
%Kinv  = getInverseChol(cholC);
%Sigma = Sigmap - (A*Sigmap)'*Kinv*(A*Sigmap);





return;



