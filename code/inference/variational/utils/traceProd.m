function tr  = traceProd( A, B )
%TRACEPROD Summary of this function goes here
%   Returns the trace of a product of two matrices trace(A*B) 

tr = sum(sum(A'.*B));


end

