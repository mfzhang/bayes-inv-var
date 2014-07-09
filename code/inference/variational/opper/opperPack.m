function theta  = opperPack( omegaq, eta )
%OPPERPACK pack parameters into a single vector theta
%   Detailed explanation goes here
% We do theta_eta = log(eta)
%
M       = length(omegaq);
Ltheta  = 2*M;
theta   = zeros(Ltheta,1);

theta(1:M)        = omegaq;
theta(M+1:Ltheta) = eta;


%% log transform
%theta(M+1:Ltheta) = log(theta(M+1:Ltheta));


return;



