function optconf = get_default_optconf()
% gets default optimization options

optconf             = optimset('GradObj','on', 'DerivativeCheck', 'off', ...
                      'LargeScale', 'off' ,'Display', 'iter');
optconf.MaxIter     = 100;  
optconf.MaxFunEvals = 1000; 
optconf.numDiff     = 0; % 0: user-supplied
                         % 1: fwd-difference, 2: central-difference

return;

