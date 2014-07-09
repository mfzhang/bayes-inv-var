function logdet  = getLogDetEigenSymmetric( V )
%GETLOGDETEIGENSYMMETRIC Summary of this function goes here
%   Detailed explanation goes here
% Gets log determinant ased on the Eigen decomposition: E V E'
%   Gets the inverse based on the Eigen decomposition: E V E'
% V: Matrix of Eigenvalues on the diagonal
logdet = sum(log(diag(V)));

end

