function [ optconf ] = opperGetOptimizerConf( optimizer, conf )
%OPPERGETOPTIMIZERCONF Summary of this function goes here
%   Detailed explanation goes here

switch optimizer,
    case 'minFunc',
        
    %% Optimization configuration
        optconf             = optimset('GradObj','off', 'DerivativeCheck', 'off', ...
                      'LargeScale', 'off' ,'Display', 'iter');
        optconf.MaxIter     = conf.MaxIter;  
        optconf.MaxFunEvals = 1000; 
        optconf.progTol     = 1e-10;
        optconf.optTol      = 1e-10;
        optconf.numDiff     = 0; % 0: user-supplied; 1: fwd-difference, 2: central-difference
        optconf.DerivativeCheck = 'off'; % SET TO OFF FOR PRODUCTION
    
    case 'nlopt',
        %localopt.algorithm      = NLOPT_LD_LBFGS;
        %localopt.verbose        = 1;        
        %optconf.local_optimizer = localopt;
        
        optconf.algorithm        = NLOPT_LD_LBFGS;
        optconf.verbose         = 1;
  
end
        
        
return;


