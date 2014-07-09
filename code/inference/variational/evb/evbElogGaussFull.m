function [elog gradMuq gradSigmaq] = evbElogGaussFull(y, fwdFunc, muq, Sigmaq, Sigmay, S)
% EVB Full Gaussian Posterior
% Computes expected log likelihood term and its gradients
% y
% fwdFunc
% muq
% Sigmaq 
% Sigmay
% S: The Number of samples to compute the MC estimates of the expectations
N = size(y,1); % Number of observation

if (length(Sigmay) == 1) % only nose variance specified
    Sigmay = Sigmay*eye(N);
end

Lq      = getChol(Sigmaq);
Ly      = getChol(Sigmay);
Lambdaq = getInverseChol(Lq); % Lambdaq = Sigmaq^-1, can pre-compute

%% Samples from posterior 
Thetas = sampleGaussChol(muq, Lq, S);

%% Evaluates Fwd Model and LogNormal at all these samples
Fval    = zeros(N,S);
LNormal = zeros(1,S); 
for s = 1 : S
    Fval(:,s)  = feval(fwdFunc, Thetas(:,s));
    LNormal(s) = logOfNormalChol(y, Fval(:,s), Ly);
end
elog = mean(LNormal);


%% Gradients wrt muq
gradMuq = zeros(size(muq));
for s = 1 : S
    gradMuq = gradMuq + (Thetas(:,s) - muq)*LNormal(s);
end
gradMuq  = gradMuq/S;
%gradMuq = solveLinChol(Sigmaq, gradMuq);
gradMuq  = Lambdaq*gradMuq;

%% Gradients wrt Sigmaq
gradSigmaq = zeros(size(Sigmaq));
for s = 1 : S
    val = Thetas(:,s) - muq;
    gradSigmaq = gradSigmaq + val*val'*LNormal(s);
end
gradSigmaq = gradSigmaq/S;
gradSigmaq = (0.5)*(Lambdaq*gradSigmaq*Lambdaq - Lambdaq*elog); 



%% Constant term
% if (length(Sigmay) == 1) % only nose variance specified
%     ldet =  - (N/2)*log(2*pi*sigmay);
% else
%     % TODO
% end
% 
% %% Samples from posterior 
% Thetas = sampleGaussChol(muq, Lq, S);
% 
% %% Evaluates Fwd Model at all these samples
% FVAL = zeros(N,S);
% for s = 1 : S
%     FVAL(:,s) = feval(fwdFunc, Thetas(:,s));
% end
% 
% %% Computes MC estimate of expectation
% lquad = 0;
% for s = 1 : S
%     lquad = lquad + (y-FVAL(:,s))'*(y-FVAL(:,s));
% end
% lquad = -(1/(2*S))*(1/sigmay)*lquad;
% 
% %% The value of the expected log likelihood
% elog = ldet + lquad;


%% Compute the gradients here


return;



