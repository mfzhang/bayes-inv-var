function [ fwdval, J, err_J ] = wrapNonlinearFwd( fwdfunc, x )
%WRAPNONLINEARFWD Summary of this function goes here
%   Detailed explanation goes here
% Wraps non-linear fwd model to also estimate the Jacobian
% - fwdfunc: Pointer to fwd model function that receives a single set of 
% parameters  x as arguments
%  -fwdval: the value(s) of the forward model evaluated at x
% J: the estimated Jacobian
% err_J: The estimated error on the Jacobian


if (nargout == 1) % Check if gradients are required
    fwdval = feval(fwdfunc, x);  % simply computes the function
else
%     nout = nargout(fwdfunc); % FIX THIS
%     if (nout > 1) % functions provides Jacobian
%         fprintf('Using Jacobian provided by fwd function ...');
%         [fwdval J] = feval(fwdfunc, x);
%     else % Estimates Jacobian numerically
        fwdval      = feval(fwdfunc, x); 
        fprintf('Computing Jacobian numerically ... ');
        [J err_J]   = jacobianest(fwdfunc, x);
        fprintf('done\n');
%     end
end
fwdval  = fwdval(:); % column vector    



return;





