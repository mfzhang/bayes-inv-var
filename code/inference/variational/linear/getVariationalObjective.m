function nelbo = getVariationalObjective(y, fobj, ...
                                 nu0, Lambda0, muq, cholLambdaq, lambday )
%GETVARIATIONALOBJECTIVE Receives functional of fwd model
%   Detailed explanation goes here
% Lambda0: Precision of prior
% nu0: Precision-adjusted mean

[fwdval J] = feval(fobj, muq);

    
D = size(cholLambdaq,1);
sigmay = 1/lambday;

if ( sum(sum(Lambda0 - diag(diag(Lambda0)))) == 0) % matrix is diagional
    Sigmap = diag(1./diag(Lambda0));
else
    cholLambdap = chol(Lambdap, 'lower');
    Sigmap      = cholLambdap'\(cholLambdap\eye(D));
end
mup    = Sigmap*nu0;

Sigmaq = cholLambdaq'\(cholLambdaq\eye(D));

nelbo = negElboLinearGaussGaussIso(y, fwdval, J, mup, Sigmap, muq, Sigmaq, sigmay);

return;


 