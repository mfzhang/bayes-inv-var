function [post_param, nlog] = learn_regbayes(y, rbconf, optconf, param)
% learns posterior parameters through regularized variational inference

%% Initialization
if (isempty(param.post))
    theta = feval(rbconf.initparamfunc, rbconf.D);
else
    theta = feval(rbconf.param2vecfunc,param.post);
end

%% optimization
fobj = @(xx) regbayes_negelbo(xx, y, rbconf, param);
[theta_new nlog] =  feval(rbconf.optfunc, fobj, theta, optconf);

%% converting optimal to standard form
post_param = feval(rbconf.vec2paramfunc, theta_new);



return;


