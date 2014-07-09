function [nlog ngrad] = regbayes_negelbo(theta, y, rbconf, param)
% wrapper for the negative of regbayes_elbo
% when posterior is a gaussuan diagonal
% TODO: put A and sigmay into the parameters as well
% As the posterior is a diagonal Gaussian we have
% L = length(theta), lenth(mu1)=L/2 and length(diagSigma1) = L/2

post_param = feval(rbconf.vec2paramfunc, theta);

if (nargout == 2)
    [log_elbo grad_elbo] =  regvarbayes_elbo( y, rbconf.klfunc,  ...
                            rbconf.elikefunc, ...
                            rbconf.regfunc, rbconf.fwdfunc, ...
                            param.prior, param.like, ...
                            post_param, param.fwd, param.reg);
    nlog  = - log_elbo;
    ngrad = - grad_elbo; 
else                       
    nlog = - regvarbayes_elbo( y, rbconf.klfunc, rbconf.elikefunc, ...
                           rbconf.regfunc, rbconf.fwdfunc, ...
                           param.prior, param.like, ...
                           post_param, param.fwd, param.reg);
end


                       
return;



 