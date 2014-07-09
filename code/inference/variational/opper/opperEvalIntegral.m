function  val  = opperEvalIntegral( f, sigma2 )
%OPPEREVALINTEGRAL Evaluate 1D integrals for Opper's method

% IT was like this before but has numerical errors
val = integral(f,-Inf,Inf);

%sigma = sqrt(sigma2);
%val = integral(f,-10*sigma,10*sigma);


return;



