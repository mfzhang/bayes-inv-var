% testVarlinearGravityGeology.m
% Tests the variational linearization approach to learn the geology in 
% a gravity problem 
% note parameters determining the geology are non-linearly transformed
% thetaGranites are the points that determine the geology (depths)
clear all; clc; close all;
randn('seed',51);

str_init = datestr(now);
diary(['results_test_varlinear_gravity_geology_', str_init, '.log']);

addpath(genpath('~/Dropbox/Matlab/gpml-matlab')); % solve_chol
addpath(genpath('~/Dropbox/Matlab/DERIVESTsuite')); % For Jacobian computations
addpath(genpath('~/Dropbox/Matlab/minFunc_2012'));
addpath(genpath('~/Dropbox/Matlab/autodiff'));

%% General settings
%RESULT_FILE     = 'tmp_50_24-Jan-2014 10:35:49.mat';
%RESULT_FILE     = 'results_test_varlinear_gravity_geology_24-Jan-2014 13:03:18.mat';
%RESULT_FILE      = 'results_test_varlinear_gravity_geology_03-Feb-2014 09:37:03.mat';
%RESULT_FILE = 'results_test_varlinear_gravity_geology_03-Feb-2014 12:00:05.mat';
%RESULT_FILE     = 'results_test_varlinear_gravity_geology_10-Feb-2014 15:07:18.mat';
RESULT_FILE = 'results_test_varlinear_gravity_geology_10-Feb-2014 17:25:41.mat';
LEARN_POSTERIOR = 1; % 0: loads from file
MAXITER         = 100; % Maximum number of iterations
TOL             = 1e-3;
nGravSens       = 100; %Approximate number of sensor locations in survey

%set(0,'DefaultFigureWindowStyle','docked');

%% Generate world
[geologySettings fixedParams thetaGranites transitions propertyFunction] = genWorld();
regionLatLon = [-32 135];
% Granites as a vector:
thetaGranites = thetaGranites{1}(:)';
thetaGranites(thetaGranites<-0.5)=-0.5;
NPARAM  = length(thetaGranites);
dimCtrl = sqrt(NPARAM);

%% Gravity survey Simulation Settings
gridSize = [64 64 64]; % Resolution of the 2-D integration

%% Generate sensor locations
%randn('seed',1); rand('seed',5);
sensorGravDepth = -0.5; % Negative for about ground. Measured in metres
edging = (gridSize-2)./gridSize;
sensorGravLocations = [rand(1,nGravSens)*geologySettings.boundaries(2)*edging(1);...
    rand(1,nGravSens)*geologySettings.boundaries(4)*edging(2); sensorGravDepth*ones(1,nGravSens)];
sensorGravLocations = [(sensorGravLocations(1,:)+geologySettings.boundaries(2)/gridSize(1));...
    sensorGravLocations(2,:)+geologySettings.boundaries(4)/gridSize(2);...
    sensorGravLocations(3,:)];

%% Settings for gravity forward model
[gravSettings] = wbt.fwd.gravmag.sensorArray3D(geologySettings.boundaries,...
    sensorGravLocations, [], gridSize, regionLatLon,2);
%theta2Transitions = @(query, ctrlPts)...
%    transitions(query, [fixedParams reshape(ctrlPts,dimCtrl,dimCtrl)]);
theta2Transitions = @(query,ctrlPts) wrapTransitions(query, ctrlPts, transitions, fixedParams, dimCtrl);

%% Grav forward model parameterised by granite vector:
theta2Grav= @(ctrlPts)wbt.fwd.gravmag.cells3D(...
     @(query)theta2Transitions(query, ctrlPts), propertyFunction, gravSettings)';
 
%%  True likelihood  with Noise Standard deviation for data generation
 data.gravFwdNoiseStd = 1e+3;
 data.gravGroundTruth = theta2Grav(thetaGranites)+randn(1,nGravSens)*data.gravFwdNoiseStd;


 %% Visualisation of true density and gravity anomaly
%  figure;
%  transitionFunction = @(query)theta2Transitions(query, thetaGranites);
%  visGravWorld(geologySettings,transitionFunction,sensorGravLocations,data.gravGroundTruth);
%  subplot(1,2,1); title('True Density'); subplot(1,2,2);  title('Gravity Observations');
h = visualizeGeology( geologySettings, sensorGravLocations, theta2Transitions, thetaGranites, data.gravGroundTruth );
set(gcf, 'Name', 'True Data');


%% Sets forward model
rbconf.fwdfunc = theta2Grav;


%% The prior in the world model is the prior used the generate the data
% but not my assumed prior for inference
% so I can use a prior here and then generate a config file 
% that parallel tempering will understand (For comparison with MCMC)
%% Canonical parameters of prior distribution: nu and precision
rbconf.D       = length(thetaGranites);
mu_p            = ones(rbconf.D,1) + 0.1*randn(rbconf.D,1); % prior is zero mean unit variance
sigma_p         = 1;
lambda_p        = 1/sigma_p;                 % isotropic gaussian
Lambda_p        = lambda_p*eye(rbconf.D);
nu_p            = Lambda_p*mu_p;
param.prior    = {nu_p; Lambda_p};

%% Here we visualize our prior mean
data_prior = theta2Grav(mu_p) + randn(1,nGravSens)*sigma_p;
%transitionFunction = @(query)theta2Transitions(query, mu_p);
%figure, visGravWorld(geologySettings,transitionFunction,sensorGravLocations,data_prior);
%subplot(1,2,1); title('Prior Density Mean'); subplot(1,2,2);  title('Prior Gravity Sample');
h = visualizeGeology( geologySettings, sensorGravLocations, theta2Transitions, mu_p, data_prior );
set(gcf, 'Name', 'Prior');



%% Parameters of likelihood
%sigmay = 1e-4*(var(data.gravGroundTruth)); % assume low nois observations
sigmay  = 1;
lambday = 1./sigmay;
param.like = {lambday};

%% empty posterior parameters force default initialization
param.post = {};

%% executes variational inference
if (LEARN_POSTERIOR)
    rbconf.maxiter = MAXITER;
    rbconf.tol     = TOL;
    y              = data.gravGroundTruth';
    [ param, muPred, sigmaPred, nelbo ] = varLinearIterativeGaussAlternate(y, param, rbconf);
%    [ param, muPred, sigmaPred, nelbo, err ] = varLinearIterativeGauss(y, param, rbconf);
else
    load (RESULT_FILE, 'muPred', 'sigmaPred', 'nelbo');
    %
    %load (RESULT_FILE, 'mu_q');
    %muPred = mu_q;
end


%% Visualize the predicted anomaly
data_post = theta2Grav(muPred);
%transitionFunction = @(query)theta2Transitions(query, muPred);
%figure, visGravWorld(geologySettings,transitionFunction,sensorGravLocations,data_post);
%subplot(1,2,1); title('Posterior Density Mean'); subplot(1,2,2);  
%title('Posterior Gravity Sample');
h = visualizeGeology( geologySettings, sensorGravLocations, theta2Transitions, muPred, data_post );
set(gcf, 'Name', 'Posterior Mean');

%% "Predicted" vs groundtruth
figure; scatter(data.gravGroundTruth, data_post); 
xlabel('Ground Truth'); ylabel('Predictions at the mean posterior');

%% actual control points values
figure; scatter(thetaGranites', muPred); 
xlabel('True Ctrl Points'); ylabel('Posterior Mean Ctrl Points');


%% Visualize the standard deviation
stdPred = sqrt(sigmaPred);
data_post = theta2Grav(stdPred);
%transitionFunction = @(query)theta2Transitions(query, stdPred);
%figure, visGravWorld(geologySettings,transitionFunction,sensorGravLocations,data_post);
%subplot(1,2,1); title('Posterior Density Std Dev'); subplot(1,2,2);  
%title('Posterior Gravity Sample');
h = visualizeGeology( geologySettings, sensorGravLocations, theta2Transitions, stdPred, data_post );
set(gcf, 'Name', 'Posterior Std');

%% Plots negative lower bound
figure, plot(1:length(nelbo), nelbo);
xlabel('Iterations'); ylabel('Negative Evidence Lower Bound');

str_end = datestr(now);
if (LEARN_POSTERIOR)
    save(['results_test_varlinear_gravity_geology_',str_end, '.mat']);
end

diary off;

