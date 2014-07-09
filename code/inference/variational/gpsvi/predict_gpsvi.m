function [mu var] = predict_gpsvi(Kmminv, covfunc, loghyper, m, S, z, xstar)
% Makes predictions with gpsvi
% z are the inducing point lcoations
% xstar are the test-pints

[Kss Kms] = feval(covfunc, loghyper, z, xstar);
Ksm       = Kms';

%% Mean of predictive distribution
mu          = Ksm*Kmminv*m;

%% variance of predictive distribution
% we can also compute full covariance at a higher cost
% diag(Ksm * kmminv * S * Kmmonv *Kms) 
var_1 =  sum(Kms.*(Kmminv*S*Kmminv*Kms),1)';
var_2 =  sum(Kms.*(Kmminv*Kms),1)';   

var = var_1 + Kss - var_2;

return;

