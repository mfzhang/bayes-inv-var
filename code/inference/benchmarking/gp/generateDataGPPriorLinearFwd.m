function  [data param] = generateDataGPPriorLinearFwd( N, PTRAIN )
%GENERATEDATAGPPIRORLINEARFWD Generate data with a GP prior and Lienar fwd
%model

FONTSIZE = 16;

   
if (nargin == 0)
    N = 100;
end


%% Settings Here
covFunc    = 'covSEard';
x          = linspace(-5,5,N)';
D          = size(x,2);
nhyper     = eval(feval(covFunc));
loghyper   = log(ones(nhyper,1));

%% Likelihood Parameters
A        =  10*rand(1,1);
fwdFunc  = @isoLinearFwdModel;
fwdParam = A;
sigma_y  = 1e-5;
lambda_y = 1/sigma_y;
Sigma_y  = sigma_y*eye(N);


%% Generate Latent Functions
f          = sampleGP(x, covFunc, loghyper, 1e-7);
plot(x,f, '.', 'MarkerSize', 12); 
xlabel('x', 'FontSize', FONTSIZE);
ylabel('f', 'FontSize', FONTSIZE);
title('Latent Function from GP prior', 'FontSize', FONTSIZE); 



%% Draws from conditional likelihhood
g  = feval(fwdFunc, f, fwdParam);
y  = sampleGauss(g, Sigma_y, 1);
figure, plot(f, y, '.', 'MarkerSize', 12); 
xlabel('Latent Function f', 'FontSize', FONTSIZE); 
ylabel('Observations y = g(f) + Noise','FontSize', FONTSIZE); 
title('Data','FontSize', FONTSIZE);

param.meanFunc     = 0; 
param.covFunc      = covFunc;
param.like.lambda  = lambda_y;
param.fwd{1}       = fwdParam;
param.fwdFunc      = fwdFunc;
param.loghyper     = loghyper;


%% Divides x into training and testing
NTRAIN    = round(N*PTRAIN);
idx       = randperm(N);
idx_train = idx(1:NTRAIN);
idx_test  = idx(NTRAIN+1:N);
data.xtrain    = x(idx_train,:);
data.ftrain    = f(idx_train,:);
data.ytrain    = y(idx_train,:);

%% Just all the data
%data.xtest     = x(idx_test,:);
%data.ftest     = f(idx_test,:);
%data.ytest     = y(idx_test,:);
data.xtest     = x;
data.ftest     = f;
data.ytest     = y;


return;



 