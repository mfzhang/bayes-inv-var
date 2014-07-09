
function energy = targetDistrib(params,cfg)

 
% Forward model
samplesGrav = cfg.modelGrav(params); 

% Prior over parameter
logprior = logpdf_prior( params, cfg );

% Unnormalised negative log posterior: 
energy = -wbt.inference.gaussianLogLikelihood(cfg.data.gravGroundTruth, ...
                             samplesGrav, cfg.data.gravFwdNoiseStd) - logprior;
                         
end


