function [ optconf ] = get_optconf_varlinear()
%GET_OPTCONF_VARLINEAR Summary of this function goes here
%   Detailed explanation goes here


optconf             = optimset('GradObj','off', 'DerivativeCheck', 'off', ...
                      'LargeScale', 'off' ,'Display', 'iter');
optconf.MaxIter     = 100;  
optconf.MaxFunEvals = 1000; 
optconf.numDiff     = 2; % 0: user-supplied
                         % 1: fwd-difference, 2: central-difference
                         
                         

end

