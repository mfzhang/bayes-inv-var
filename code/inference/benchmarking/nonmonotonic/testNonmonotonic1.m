clear all; clc;

%% Prior parameters
M      = 5000;
mu_p    = zeros(M,1);
sigma_p = 1;
Sigma_p = sigma_p*eye(M);


%% Draws one sample from N-dimensional prior
theta   = sampleGauss(mu_p,Sigma_p,1); 

%% Non-monotonic fwd model
%fwdfunc      = @(xx) xx.^2;
fwdfunc       = @(xx) xx.*(0.5-xx).*(1-xx);
f             = feval(fwdfunc, theta);
figure,plot(theta,f, '.');

%% Conditional likelihood
y = f;

EDGES = -1:0.001:1;
N = histc(y, EDGES );
figure, bar(EDGES, N, 'histc');
