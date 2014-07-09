function muY  = opperPredictObservable(pred, fwdFunc, ...
                                       fwdParam, conf )
%OPPERPREDICTOBSERVABLE Predicts observable y = g(f)
%   Detailed explanation goes here
% pred : the preditive distribution
% fwdFunc: the fwd function g
% fwdParam: The parameters of the fwd Function
% conf.S: The number of sammples to use in case of MC estiamte
% conf.elogMethod: The integration method to use

%% Read first configuration
S          = conf.S;
fwdFunc    = @(xx) fwdFunc(xx, fwdParam{:});
elogMethod = conf.elogMethod;


switch elogMethod,
        case 'quad',
            fnElog= @opperPredGaussQuad;
        case 'mc', % TO DO
            fnElog= @opperPredGaussMC;
     otherwise,
        disp('Invalid ELOG Method');
end

muY = feval(fnElog, fwdFunc, pred.mu, pred.Sigma, S);



return;


%% Using Quaf
function muY = opperPredGaussQuad(fwdFunc, muStar, SigmaStar, S)
nTest = size(muStar,1);
muY   = zeros(nTest,1);
for i = 1 : nTest
    h = @(ff)  fPredIntegrand(ff, fwdFunc, muStar(i), SigmaStar(i,i));   % integrand
    muY(i) =  opperEvalIntegral(h, SigmaStar(i,i)) ; %integral(h,-Inf,Inf);
end


return;


%% MC Estimate
function muY = opperPredGaussMC(fwdFunc, muStar, SigmaStar, S)
nTest = size(muStar,1);
diagSigmaStar = diag(SigmaStar);
Thetas = sampleGaussDiag(muStar, diagSigmaStar, S);
fwdVal = zeros(nTest,S);
for s = 1 : S
   theta             =  Thetas(:,s);
   fwdVal(:,s)        = feval(fwdFunc, theta);  
end
muY = mean(fwdVal,2);

return;




%% Integrand for computing the predictions
function val = fPredIntegrand(theta, fwdFunc, muStar, sigmaStar)
fwdVal      = feval(fwdFunc, theta);
gaussVal    = valOfGaussUnivariate(theta, muStar, sigmaStar);
val         = fwdVal.*gaussVal;

return;