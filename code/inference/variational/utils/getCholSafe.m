function  L  = getCholSafe( Sigma, maxTries, minJitter )
%GETCHOLSAFE Summary of this function goes here
%   Detailed explanation goes here

% IT WAS LIKE THIS BEFORE
% if (nargin == 1)
%     jitter = 1e-7; % TODO: Get enough jitter automatically
% end
% 
% Sigma = Sigma + jitter*eye(size(Sigma));
% L     = getChol(Sigma);
% 
% return;

if (nargin < 3)
    minJitter = 1e-6;
end

if nargin < 2
  maxTries = 10;
end
jitter = 0;
for i = 1:maxTries
  try
    % Try --- need to check Sigma is positive definite
    if jitter == 0;
      jitter = abs(mean(diag(Sigma)))*minJitter;
      L = getChol(Sigma);
      break
    else
      if nargout < 2
        warning(['Matrix is not positive definite in jitChol, adding ' num2str(jitter) ' jitter.'])
      end
      L = getChol(real(Sigma+jitter*eye(size(Sigma, 1))));
      break
    end
  catch
    % Was the error due to not positive definite?
    nonPosDef = 0;
    verString = version;
    if str2double(verString(1:3)) > 6.1
      [void, errid] = lasterr;
      if strcmp(errid, 'MATLAB:posdef')
        nonPosDef = 1;
      end
    else
      errMsg = lasterr;
      if findstr(errMsg, 'positive definite')
        nonPosDef = 1;
      end
    end
  end
  if nonPosDef
    jitter = jitter*10;
    if i==maxTries
      error(['Matrix is non positive definite tried ' num2str(i) ...
             ' times adding jitter, but failed with jitter ' ...
             'of ' num2str(jitter) '. Increase max tries'])
    end
  else
    error(lasterr)
  end
end



return;

