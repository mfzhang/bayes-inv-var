function [ omegaq, eta ] = opperUnpack( theta )
%UNPACKOPPER unpack parameters 
%   Detailed explanation goes here
% we do: eta = exp(theta_eta)

Ltheta = length(theta);
M      = Ltheta/2; % M parameters for the posterior parameters
         % and M parameters for the "free" diagonal parameters for the covariance

omegaq     = theta(1:M);
eta        = theta(M+1:Ltheta); 

%% exp transform        
%eta  = exp(eta);         


return;

