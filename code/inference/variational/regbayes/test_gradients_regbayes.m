%% checks the derivatives
function test_gradients_regbayes()
path(path(), genpath('~/Dropbox/Matlab/DERIVESTsuite'));

N          = 200;
D          = 2;
data_fname = 'toydata.mat';
generate_toy_data(200, D, data_fname);
check_gradients(N, D, data_fname);



%% 
function check_gradients(n, D, data_fname)
y = []; mu0 = []; sigma0 = []; sigmay = []; A = []; 
load(data_fname, 'y', 'mu0', 'sigma0', 'sigmay', 'A');
y = y(1:n);
A = A(1:n,:);


%% algorithm configuration
rbconf.klfunc         =  @kl_gauss_iso_diag;
rbconf.elikefunc      =  @eloglike_gauss_iso_diag;
rbconf.fwdfunc        = []; % forward model
rbconf.regfunc        = []; % Regularizer
rbconf.param2vecfunc  =  @param2vec_gauss_diag;
rbconf.vec2paramfunc  =  @vec2param_gauss_diag;
rbconf.initparamfunc  =  @initparam_gauss_diag;
rbconf.optfunc        =  @minFunc; % {'minFunc', 'fminunc'}
rbconf.D              =   D; % Number of parameters

%% parameters of distributions
param.prior = {mu0, sigma0};
param.like  =  {A, sigmay};
param.fwd   = [];

myfunc = @(xx) regbayes_negelbo(xx, y, rbconf, param);

R = 10;
L = 2*D;
delta = zeros(L,R);
for r = 1 : 10
    theta  = rand(L,1);
    [grad_diff, errest, finaldelta] = gradest(myfunc, theta);
    [nlog grad_func] = myfunc(theta);
    delta(:,r) = grad_diff' - grad_func;
    
    delta(:,r)
end
hist(delta(:))


% derivativeCheck(@regbayes_negelbo, xx, 1, type, y, rbconf, param );

