function [ fwdval, J, err_J ] = wrapNonlinearFwdBatch( fwdfunc, x, idx )
% evaluates fwd model and returns  subset of observations 
% given by index idx


ptrfunc = @(xx) wrapFwdBatch(fwdfunc, xx, idx);

if (nargout == 1) % Check if gradients are required
    fwdval = feval(ptrfunc, x);  % simply computes the function
else
%     nout = nargout(fwdfunc); % FIX THIS
%     if (nout > 1) % functions provides Jacobian
%         fprintf('Using Jacobian provided by fwd function ...');
%         [fwdval J] = feval(fwdfunc, x);
%     else % Estimates Jacobian numerically
        fwdval      = feval(ptrfunc, x); 
        fprintf('Computing Jacobian numerically ... ');
        [J err_J]   = jacobianest(ptrfunc, x);
        fprintf('done\n');
%     end
end
fwdval  = fwdval(:); % column vector    



return;


function  fwdval  = wrapFwdBatch( fwdfunc, x, idx )
%WRAPFWDBATCH evaluates fwd model and returns  subset of observations 
% given by index idx

fwdval = feval(fwdfunc, x);
fwdval = fwdval(idx);


return;







