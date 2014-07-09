%% Initial settings for toy test
function y = generate_toy_data(N, D, data_fname)
mu0        = zeros(D,1); % prior
sigma0     =    0.01;    % parameters
%mu1        = rand(D,1);       % approximate posterior
%diagSigma1 = 0.01*rand(D,1);  % parameters
sigmay     = 0.01;       % likelihood
A          = rand(N, D); % parameters

%% Draws from joint prior-likelihood model
theta      = mu0 + sqrt(sigma0)*randn(D,1);
y          = A*theta + sqrt(sigmay)*randn(N,1);

save(data_fname, 'y', 'mu0', 'sigma0', 'sigmay', 'A', 'theta');

return;


