function nelbo = getVariationalObjectiveFixed(y, fwdval, J, ...
                                 nu0, Lambda0, muq, cholLambdaq, lambday )
%GETVARIATIONALOBJECTIVEFIXED Summary of this function goes here
%   Detailed explanation goes here
%% getVariationalObjectiveFixed
% gets variationa objective for fixed fwdval and Jacobian

% Lambda0: Precision of prior
% nu0: Precision-adjusted mean
    
D = size(cholLambdaq,1);
sigmay = 1/lambday;

if ( isDiag(Lambda0)) % matrix is diagional
    Sigmap = diag(1./diag(Lambda0));
else
    cholLambdap = chol(Lambdap, 'lower');
    Sigmap      = cholLambdap'\(cholLambdap\eye(D));
end
mup    = Sigmap*nu0;

Sigmaq = cholLambdaq'\(cholLambdaq\eye(D));

nelbo = negElboLinearGaussGaussIso(y, fwdval, J, mup, Sigmap, muq, Sigmaq, sigmay);

return;




