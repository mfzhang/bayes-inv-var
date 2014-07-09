function [y, J, errJ] = simple_grav_fwd(rho, G)
% simple deterministic gravity fwd model
% 07/01/2014: Changed order of the imput so that the parameters are 
% at the end
y = G*rho;

if (nargout >= 2) % evaluates the Jacobian at rho 
    J = G;
    errJ = 0;
end
    
return;
