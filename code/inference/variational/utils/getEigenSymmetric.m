function [ E V ] = getEigenSymmetric( K, minNoise )
%GETEIGENSYMMETRIC Summary of this function goes here
%   Gets egeindecomposition of a symmetric matrix K
D = size(K,1);

if (nargin==1)
    minNoise = 0;
end

[E V] = eig(K + minNoise * eye(D,D)); 

%% Gets rid off possible small complex and negatives
V(V<0) = 0;  
V      = real(V); 





end

