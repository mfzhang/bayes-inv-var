function  [elog gradMuq, gradSigmaq]  = opperElogGaussFullMC(y, fwdFunc, ...
                                         muq, Sigmaq, diagSigmay, S )
%OPPERELOGGAUSSFULL Expected log likelihood for Opper's model
%   Detailed explanation goes here
% diagSigmay: either a number or a vector with the diagonal values
% fwd functions acts on a per single point basis but can also 
% evaluate N (independent) points at a time

%% DELETE ME?
rng('default'); 

N = size(y,1);

if (length(diagSigmay) == 1) % only noise variance specified
    diagSigmay = diagSigmay*ones(N,1);
end

%% Gets the diagonal of the covariance as we only need the marginals
diagSigmaq    = diag(Sigmaq); % Vector representation of diagonal of Sigmaq


%% Samples from current marginal posterior 
Thetas = sampleGaussDiag(muq, diagSigmaq, S);

%% Evaluates Fwd Model and LogNormal at all these samples
Fval    = zeros(N,S);
LNormal = zeros(N,S); 
for s = 1 : S
    Fval(:,s)  = feval(fwdFunc, Thetas(:,s));
    LNormal(:,s) = logOfGaussUnivariate(y, Fval(:,s), diagSigmay); 
end
velog = mean(LNormal,2); % I will use this one later
elog  = sum(velog); 

%% Computes the gradients if required
if (nargout > 1)
    %% Gradients wrt muq
    gradMuq = zeros(size(muq));
    for s = 1 : S
     gradMuq = gradMuq + (Thetas(:,s) - muq).*LNormal(:,s);
    end
    gradMuq  = gradMuq/S;
    gradMuq  = gradMuq./diagSigmaq;


    %% Gradients wrt Sigmaq is diagonal (we only need a vector to represent it)
    gradSigmaq = zeros(size(Sigmaq,1),1);
    for s = 1 : S
        val = Thetas(:,s) - muq;
        gradSigmaq = gradSigmaq + (val.^2).*LNormal(:,s);
    end
    gradSigmaq = gradSigmaq/S;
    gradSigmaq = (0.5)*(gradSigmaq - diagSigmaq.*velog)./(diagSigmaq.^2); 
end


%% DELETE ME
% elog    = muq'*muq;
% gradMuq = 2*muq;
% gradSigmaq = zeros(size(gradSigmaq));

return;


