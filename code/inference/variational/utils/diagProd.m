function v = diagProd(A, B)
% computes v = diag(A*B)
% A is a nxm matrix and B is a mxn matrix 
% so that the resulting matric C = A*B is squared

v = sum(A'.*B,1)';


return;

