function x = solveLinChol( L, y )
%SOLVELINCHOL Solves linear sytem Ax = y using chol decom of A
%   Detailed explanation goes here
% L: Lower triangular chol
x = solve_chol(L',y);

return;

