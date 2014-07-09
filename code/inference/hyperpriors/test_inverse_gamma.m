% tests of inverse-gamma priors for noise variance parameters
clear all; clc;

%% General settings
v_alpha = [0.01 0.1 0.1 1 2 1 2] ; % shape parameter 
v_beta  = [0.01  1  2   1 1 2 2]; % rate parameter
N     = 1000;


%% x ~ Gamma(alpha,1/beta) -> y ~ InvGamma(alpha, beta)
L = length(v_alpha);
Y = zeros(N,L)
for i = 1 : L
    x = gamrnd(v_alpha(i), 1/v_beta(i), N, 1);
    Y(:,i) = 1./x;
end
cdfplot(Y);

