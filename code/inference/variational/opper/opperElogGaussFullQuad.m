function [elog gradMuq, gradSigmaq] = opperElogGaussFullQuad( y, fwdFunc, ...
                                         muq, Sigmaq, diagSigmay, S )
%OPPERELOGGAUSSFULLQUAD Elog for Opper's method using quadrature
%   Detailed explanation goes here
% S argument is irrelevant here but we need it for consistency with 
% opperElogGaussFullMC
%
N = size(y,1);

%% Cmputes ELOG using quadrature
velog = zeros(N,1);
for i = 1 : N
    f = @(xx) fElogIntegrand(xx, fwdFunc, y(i), muq(i), Sigmaq(i,i), diagSigmay(i));
    velog(i) = opperEvalIntegral( f, Sigmaq(i,i) ); % integral(f,-Inf,Inf);
end
elog = sum(velog);

if (nargout > 1)
    gradMuq    = zeros(size(muq));
    gradSigmaq = zeros(size(Sigmaq,1),1);    
    for i = 1 : N
        %% Gradients wrt Muq
        f = @(xx) fGradMuIntegrand(xx, fwdFunc, y(i), muq(i), Sigmaq(i,i), diagSigmay(i));
        valIntegral = opperEvalIntegral( f, Sigmaq(i,i) ); %integral(f,-Inf,Inf);
        gradMuq(i)  = valIntegral/Sigmaq(i,i); 
        
        %% Gradients wrt diag(Sigmaq)
        f = @(xx) fGradSigmaIntegrand(xx, fwdFunc, y(i), muq(i), Sigmaq(i,i), diagSigmay(i));
        valIntegral   = opperEvalIntegral( f, Sigmaq(i,i) ); %integral(f,-Inf,Inf);
        gradSigmaq(i) = (0.5)*( valIntegral - Sigmaq(i,i)*velog(i) )/(Sigmaq(i,i))^2;
    end
end



return;




%% Integrand for ELOG computation
function val = fElogIntegrand(theta, fwdFunc, y, mq, sigmaq, sigmay)

fwdVal      = feval(fwdFunc, theta);
logGaussVal = logOfGaussUnivariate(y, fwdVal, sigmay);
GaussVal    = valOfGaussUnivariate(theta, mq, sigmaq);

val = logGaussVal.*GaussVal;

return;


%% Integrand useful for comptuation of gradient wrt mu
function val = fGradMuIntegrand(theta, fwdFunc, y, mq, sigmaq, sigmay)
fwdVal      = feval(fwdFunc, theta);
logGaussVal = logOfGaussUnivariate(y, fwdVal, sigmay);
GaussVal    = valOfGaussUnivariate(theta, mq, sigmaq);

val         = (theta - mq).*logGaussVal.*GaussVal;

return;


%% Integrand useful for the computation of gradient wrt Sigma
function val = fGradSigmaIntegrand(theta, fwdFunc, y, mq, sigmaq, sigmay)
fwdVal      = feval(fwdFunc, theta);
logGaussVal = logOfGaussUnivariate(y, fwdVal, sigmay);
GaussVal    = valOfGaussUnivariate(theta, mq, sigmaq);


val         = ((theta - mq).^2).*logGaussVal.*GaussVal;

return;

