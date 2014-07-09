function  Kinv  = getInverseChol( L )
%GETINVERSECHOL Get inverse of matrix with Cholesky Decompositon
%   L: The lower triangular representation of the matrix

I    = eye(size(L));
Kinv = L'\(L\I);

% Using hte mex file doesn't work as solve_chol solves for each of the
% columns indepedently
% Kinv = solve_chol(L,solve_chol(L',I)); % solve_chol receives upper triangular


end

