function  C  = AdiagB( A, B )
% A*diag(B)

if size(B,2) == 1 % diag(A) is given
    C = A.*repmat(B',size(A,1),1);
else
    C = A.*repmat(diag(B)',size(A,1),1);
end

return;

