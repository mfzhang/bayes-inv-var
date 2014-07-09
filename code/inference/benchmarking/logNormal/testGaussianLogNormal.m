% Gaussian prior, log normal fwd model
% --> Exact Gaussian posterior
clear all; clc; 


R = 10;
muPred    = zeros(R,1);
sigmaPred = zeros(R,1);
muTrue    = zeros(R,1);
sigmaTrue = zeros(R,1);
for i = 1 : R
    %% Get data and parameters of generation process
    [y param]  = generateDataLogNormal(1000);
    param.post = {};


    % % %% sets configuration
    conf.maxiter = 2; % IN PRINCIPLE, WE DON'T NEED MORE THAN 2 GLOBAL ITERATIONS
    conf.tol     = 1e-3;
    conf.fwdfunc = @(xx) logNormalFwdModel(xx, param.fwd{:});
    % 
    % %% Runs variational
    [ param, muPred(i), sigmaPred(i)] = varLinearIterativeGaussAlternate(y, param, conf);


    % %% True posterior
    [muTrue(i), sigmaTrue(i)] = getExactPosteriorLogNormal(param,y);
end

%% Prints error results
figure;
scatter(muTrue, muPred); xlabel('True Mu'); ylabel('Pred Mu');
figure;
scatter(sigmaTrue, sigmaPred); xlabel('True var'); ylabel('Pred var');
set(gca, 'yscale', 'log', 'xscale', 'log');

