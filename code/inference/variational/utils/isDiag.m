function  val  = isDiag( K )
%ISDIAG Test if matrix K is diagonal

val =  sum(sum(K - diag(diag(K)))) == 0;



return;

