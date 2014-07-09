% test_regbayes1.m
% test regbayes with toy data
function test_regbayes1()

%clear all; clc; close all;

%% paths
path(path(), genpath('~/Dropbox/Matlab/minFunc_2012'));
path(path(), genpath('~/Dropbox/Matlab/utils'));

%rand('seed', 12);
%randn('seed', 24);

D     = 2;  % dimensionality of parameter space
MAX_N = 2000;

data_fname = 'toydata.mat';
generate_toy_data(MAX_N, D, data_fname);


%% runs inference for several data sizes
%N = 10; % Number of observations
vec_n = [100 200 500 1000 2000];
L = length(vec_n);
err_mean = zeros(1,L);
for i = 1 : L
    n = vec_n(i);
    [post_param true_param ] = test_regbayes_single(n, D, data_fname);
    err_mean(i) = get_errors(post_param, true_param);
end
plot(vec_n, err_mean, '-o');


[post_param{1} true_param{1}]




%% simple Gaussian-Gaussian test
function [post_param, true_param] = test_regbayes_single(n, D, data_fname)

y = []; mu0 = []; sigma0 = []; sigmay = []; A = []; 
load(data_fname, 'y', 'mu0', 'sigma0', 'sigmay', 'A');
y = y(1:n);
A = A(1:n,:);



% %%
%ptr_kl     = @kl_gauss_iso_diag;
%ptr_elike  = @eloglike_gauss_iso_diag;


% %% basic function
% kl = kl_gauss_iso_diag(mu0, sigma0, mu1, diagSigma1);
% elog = eloglike_gauss_iso_diag(y, A, sigmay, mu1, diagSigma1);



%prior_param = {mu0, sigma0};
%like_param  = {A, sigmay};
%post_param  = {mu1, diagSigma1};
%kl_param = [prior_param post_param];
%kl = feval(ptr_kl, kl_param{:});
%elog_param  = [like_param, post_param];
%elog        = feval(ptr_elike, y, elog_param{:});
%elbo = regvarbayes_elbo(y, ptr_kl, ptr_elike, [], [], ...
%                        prior_param, like_param, post_param, [] );
%theta  = [mu1(:) ; diagSigma1(:)];
%nlog = regbayes_negelbo_gauss_diag(theta, y, ptr_kl, ptr_elike, [], [], ...
%                        prior_param, like_param, [] ) ;


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


%% optimizer configuration
optconf = get_default_optconf();


%% parameters of distributions
param.prior = {mu0, sigma0};
param.like  =  {A, sigmay};
param.fwd   = [];


%% runs variational inference
[post_param, nlog] = learn_regbayes(y, rbconf, optconf, param);
mu_pred            = post_param{1};
sigma_pred         = post_param{2};

                    
%% true posterior
Sigma0     = sigma0*eye(D);
Sigmay     = sigmay*eye(n);
C          = A*Sigma0*A' + Sigmay;
Lc         = chol_safe(C, 1e-7);
invC       = invChol(Lc);
mu_true    = mu0 + (A*Sigma0)'*invC*(y - A*mu0);
Sigma_true = Sigma0 - (A*Sigma0)'*invC*(A*Sigma0);
true_param{1} = mu_true;
true_param{2} = Sigma_true;


return;

%% compute errors on parameter estimates
function err_mean = get_errors(post_param, true_param)

%% Compare predicted vs true means
err_mean = norm(true_param{1}-post_param{1},2)/norm(true_param{1},2);



return;















