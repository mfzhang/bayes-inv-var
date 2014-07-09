function theta = initparam_gauss_diag(D)

Mu        = rand(D,1);
diagSigma = 0.01*rand(D,1);

param = {Mu, diagSigma};

theta = param2vec_gauss_diag(param); 

return;