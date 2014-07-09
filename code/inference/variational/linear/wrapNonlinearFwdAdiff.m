function [ fwdval, J ] = wrapNonlinearFwdAdiff( fwdfunc,  x )
%WRAPNONLINEARFWDADIFF Wraps nonlinear fwd model with automatic differentiation
%   Detailed explanation goes here
% 
adObj  = adiff(x);
fwdObj = feval(fwdfunc, adObj);

fwdval = adiffget(fwdObj,'value');

if (nargout > 1)
    fprintf('Computing Jacobian with automatic differentiation... ');
    J = adiffget(fwdObj,'derivative');
    fprintf('done\n');

end

return;



