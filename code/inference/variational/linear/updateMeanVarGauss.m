function  nu_q  = updateMeanVarGauss( mu_q, nu_p, fwdval, J, y, lambday, ...
                                    cholL_q, niter, tol)
%UPDATEMEANVARGAUSS Updates the mean iteratively with closed-form updates
%   Detailed explanation goes here

fprintf(' *** Updating Mean ***\n');

D = length(mu_q);
err = zeros(niter,1);
i = 0;
mse = inf;
while ( (i < niter) && (mse > tol) )
 muOld = mu_q;
 nu_q     = nu_p + J'*lambday*(y - fwdval + J*mu_q);
 mu_q     = solve_chol(cholL_q',nu_q);       
 mse      = sum((mu_q - muOld).^2)/D;
 i = i + 1;
 err(i)  = mse;
 
 fprintf('Mean Update Iteration: %d,', i);
 fprintf('[error=%.6f]', err(i));
 fprintf('\n');
end


return;

