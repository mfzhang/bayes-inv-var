function [ omega, eta ] = opperInitParametersReparam(paramPost, D )
%OPPERINITPARAMETERSREPARAM Initializaes parameters using the
%reparamtereization directly

if (isempty(paramPost))
    omega  = zeros(D,1);    
    eta    = 1./1e-3*ones(D,1);           
else % simply assigns whatever is given
    omega  = paramPost.omega; 
    eta    = paramPost.eta;
end




return;
