function logdet = getLogDetChol( L )
%GETLOGDETCHOL Get log determinant of matrix based on Cholesky
%decomposition
%   L: The lower triangular part of the decomposition

logdet = 2*sum(log(diag(L)));

end

