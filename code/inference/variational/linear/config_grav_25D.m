function [cfg] = config_grav_25D(params)
%CONFIG configure negative feedback oscillatior and parallel tempering options
%
%   [cfg] = config()
%
% Initialize default configuration options 
cfg = init_config_defaults();
cfg.max_init_steps = 100;
% Initialise paramaters
if nargin>0
    cfg.params_init = params;
end
cfg.adapt_beta = 1; % set to 1 to adapt beta, 0 to clamp beta

%% Parameter definitions (see BNGL model for parameter descriptions) [REQUIRED]
% Available prior distributions:
%   'point'           args: 'value'
%   'uniform'         args: 'min', 'max'
%   'normal',         args: 'mu', 'sigma'
%   'lognormal',      args: 'mu', 'sigma'
%   'laplace',        args: 'mu', 'b'
%   'boundedlaplace', args: 'mu', 'b', 'min', 'max'
%   'beta',           args: 'alpha', 'beta'
%   'exponential',    args: 'mu'
%scale = 0.5;  % global weight parameter for prior strength (or "width"). scale>0.
cfg.param_defs = { ...
  struct('name','X_Param01',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param02',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param03',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param04',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param05',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param06',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param07',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param08',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param09',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param10',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param11',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param12',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param13',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param14',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param15',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param16',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param17',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param18',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param19',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param20',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param21',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param22',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param23',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param24',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  struct('name','X_Param25',  'prior','uniform', 'min', -1, 'max',2.5, 'units','(none)'  ), ...
  };
% cfg.param_defs = { ...
%   struct('name','X_Param01',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param02',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param03',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param04',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param05',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param06',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param07',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param08',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param09',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param10',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param11',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param12',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param13',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param14',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param15',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param16',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param17',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param18',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param19',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param20',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param21',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param22',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param23',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param24',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   struct('name','X_Param25',  'prior','normal', 'mu', 0, 'sigma',1, 'units','(none)'  ), ...
%   };
% initialize parameter distributions
cfg = init_parameter_defs( cfg.param_defs, cfg );


%% Observable definitions [REQUIRED]
%  Set display=0 to disable observable plot during visualization. Set
%  minplot/maxplot to see y-axis bounds in visualization scripts.
%  The 'norm' field should be the name of the observable that is used to
%  normalize the observable (e.g. YP is normalized by YT, the total quantity
%  of Y). Leave this field empty if normalization is not desired.
cfg.obsv_defs = { ...
  struct('name','XT_Param', 'units','fraction','display',1, 'minplot',-4, 'maxplot',4.0), ...
  };
% initialize observable structs
cfg = init_observable_defs( cfg.obsv_defs, cfg );


% parameters for visualizing trajectories [REQUIRED only for visualization]
cfg.time_units = 's';                               % time units
cfg.sim_tstart = 0;                                 % simulation start time, s
cfg.sim_tstop  = 63;                                % simulation stop time, s
cfg.sim_dt     = 0.5;                               % time step for trajectory, s

% penalty for long integration times [OPTIONAL: useful for avoiding stiff parameter regions]
cfg.timepenalty = 0;


%% Parallel temperating options
% Defaults are usually ok. Things you may want to change: jobname, nchains,
% parallel, nswaps, adapt_last, energy_init_max, relstep_init.
% See core/init_config_defaults.m for a complete list of config options.
cfg.jobname = 'ptemp_grav_25D';                % job name, for file input/output
cfg.parallel = 0;                      % parallel? true/false
cfg.maxlabs  = 6;                      % maximum number of labs for parallel computing
cfg.nchains  = 10;                      % number of chains
cfg.nswaps = 250100;                    % number of chain swaps (=number of saved samples!)
cfg.nsteps = 10;    %25                   % number of steps between chain swaps
cfg.display_status_interval = 10;      % How often to display info
cfg.save_progress_interval = 300;     % How often to save progress 
cfg.adapt_last = 40000;                 % last adaption step
cfg.energy_init_max = 160;             % maximum allowed energy for initialization
cfg.beta_init = 0.7;                 % beta initialization parameter
cfg.relstep_init = 0.0125;               % relstep initialization parameter


% %% Load experimental data [REQUIRED].
load data; 
cfg.data = data;
load modelGrav
% load modelMag
cfg.modelGrav = modelGrav;
% cfg.modelMag = modelMag;

%% load default functions
cfg = setup_default_functions( cfg );

% Energy function [REQUIRED!]
%   prototype: [energy] = @(params)
%
cfg.energy_fcn = @(params) targetDistribGrav(params, cfg);


% %% setup custom function handles
% 
% % Sample parameter prior function [optional]:
% %   prototype: [params] = @()
% %
% %cfg.sample_prior_fcn = @() ( <insert custom function here> );
% 
% 
% % Log-prior pdf function [optional]:
% %   prototype: [logpdf] = @(params)
% %
% %cfg.logpdf_prior_fcn = @(params) logpdf_prior( <insert custom function here> );
% 
% 
% % Equilibration protocol function [optional]:
% %   prototype: [err,sim,obsv] = @(params)
% %
% % NOTE: Equilibration is performed once per parameter set. If equilibration
% % is dependent on the protocol, include equilibration in protocol function.
% %
% %cfg.equilibrate_fcn = @(params) ( <insert custom function here> );
% 
% 
% % Define a simulation protocol for each experiment [REQUIRED!]
% %
% % pass extra options to the protocol fcn here:
% args = struct( ...
%     'param_map',    cfg.param_map,    ...% useful for finding params by name
%     'obsv_to_norm', cfg.obsv_to_norm, ...% required by "norm_obsv"
%     'obsv_norm_by', cfg.obsv_norm_by  ...% required by "norm_obsv"
% );
% % experiment protocols here:
% for d = 1 : cfg.nexpt
%     % protocol specific parameters:
%     S = cfg.data{d}.S;  % input signal to NFO
%     % Experimental protocol function:
%     %   prototype: [err,sim,obsv] = @(t,init,params)
%     cfg.data{d}.protocol_fcn = @(t,init,params) simulate(t,init,params,S,args);
% end






