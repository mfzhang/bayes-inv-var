function [ post_param nlog ] = learn_varlinear( y, rbconf, optconf, param )
%LEARN_VARLINEAR Summary of this function goes here
%   Detailed explanation goes here


%% Initialization
if (isempty(param.post))
    theta = feval(rbconf.initparamfunc, rbconf.D);
else
    theta = feval(rbconf.param2vecfunc,param.post);
end

%% optimization
fobj = @(xx) varlinear_negelbo(xx, y, rbconf, param);
[theta_new nlog] =  feval(rbconf.optfunc, fobj, theta, optconf);

%% converting optimal to standard form
post_param = feval(rbconf.vec2paramfunc, theta_new);



return;




end

