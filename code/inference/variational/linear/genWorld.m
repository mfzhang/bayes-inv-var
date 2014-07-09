function [geologySettings thetaFixed thetaGranites transitions propertyFunction] = genWorld()

%% Geology Model (use the standard demo)
[geologySettings theta transitions] = wbt.prior.setupExampleWorld();
thetaFixed = theta(1:end-1);
thetaGranites = theta(end);


% (Constant) Property Model
den = [0.0,0.0, 0.0, 1.0, 0.0];
sus = [0.0, 0.0, 0.0, 0.0, 1.0e-5];
cond = [1.0, 1.0, 1.0];
res = [1.0, 1.0, 1.0];
vel = [4500, 5000, 6000, 4000, 5300];
propertyParams = wbt.prior.ConstantPropertyParams(den, sus, cond,...
                                                  res, vel);
propertyFunction = @(rockTypes, propertyName)...
    (wbt.prior.constantPropertyModel(...
        rockTypes, propertyName, propertyParams));

    
    