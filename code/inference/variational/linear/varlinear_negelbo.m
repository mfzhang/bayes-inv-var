function nelbo  = varlinear_negelbo( theta, y, rbconf, param )
%VARLINAR_NEGELBO Summary of this function goes here
%   Detailed explanation goes here
% neg elbo for variational linear

post_param = feval(rbconf.vec2paramfunc, theta);

kl_param    = [param.prior, post_param];


kl      = feval(rbconf.klfunc, kl_param{:});
elog    = feval(rbconf.elikefunc, y, param.like{:}, rbconf.fwdfunc, ...
                    param.fwd{:}, post_param{:} );

elbo  = elog - kl;

nelbo = - elbo;
    


end


