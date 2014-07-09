function  L  = getChol( Sigma )
%GETCHOL Summary of this function goes here
%   Detailed explanation goes here

if ( isDiag(Sigma) )
    L = sqrt(Sigma);
else
    L = chol(Sigma, 'lower'); % 
end



return;
