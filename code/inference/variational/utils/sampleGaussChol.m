function  X  = sampleGaussChol( mu, L, N )
%SAMPLEGAUSSCHOL Samples from Gaussian for a given Chol decomp of Sigma
% L: The lower cholesky of Sigma
% N the numer of samples

D = length(mu);
z =  randn(D,N); % D N-dimensional independent Gaussian vectors
X = repmat(mu,1,N) + L*z; % Correlated Gaussian vectors 


return;



