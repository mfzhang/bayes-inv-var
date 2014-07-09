function  logn  = logOfNormalChol( y, m, L )
% Computes the log of a normal distribution evaluated at y
%   Detailed explanation goes here
% L: the lower cholesky decomposition of Sigma

D     = length(m);
ldet  = getLogDetChol(L);
val   = y - m;
lquad = val'*solveLinChol(L,val);

logn = -0.5*( D*log(2*pi) + ldet + lquad );


return;

