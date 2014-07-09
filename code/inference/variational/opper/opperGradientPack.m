function [ gradTheta ] = opperGradientPack( gradOmegaq, gradEta, eta )
%OPPERGRADIENTPACK Summary of this function goes here
%   Detailed explanation goes here
M                  = length(gradOmegaq);
Ltheta             = 2*M;
gradTheta          = zeros(Ltheta,1);


gradTheta(1:M)        = gradOmegaq;
gradTheta(M+1:Ltheta) = gradEta; % Accounting for log transform

%% from exp transfrom
%gradTheta(M+1:Ltheta) = gradTheta(M+1:Ltheta).*eta; 


return;
