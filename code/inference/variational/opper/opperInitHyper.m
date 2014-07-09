function  loghyper  = opperInitHyper( covFunc, x, y )
%OPPERINITHYPER Initializes hyperparameters for Opper&A
%   Detailed explanation goes here
D      = size(x,2);
nhyper = eval(feval(covFunc));
loghyper = zeros(nhyper,1);


%loghyper(1:nhyper-1) = log(1./var(x,0, 1));
%loghyper(nhyper)     = log(std(y)); % Signal variance

log(ones(nhyper,1));


return;

