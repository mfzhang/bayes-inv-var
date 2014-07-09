function [ x, fval ] = optimizeObjective(optimizer, fobj, x0, optconf )
%OPTIMIZEOBJECTIVE Summary of this function goes here
%   Detailed explanation goes here

switch optimizer,
    case 'minFunc',
        [x, fval] = minFunc(fobj, x0, optconf);    
    case 'nlopt',
        %optconf.local_optimizer.min_objective  = fobj;
        optconf.min_objective   = fobj;
        [x, fval, retcode] = nlopt_optimize(optconf,x0);
        optout.nlog     = nlog;
        optout.exitflag = retcode;        
end        






