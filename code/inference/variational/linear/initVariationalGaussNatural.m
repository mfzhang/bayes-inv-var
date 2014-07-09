function [nu Lambda cholLambda mu ] = initVariationalGaussNatural( D, paramPost )
%INITVARIATIONALGAUSSNATURAL Initializes variational approach with Gaussian posterior
%  Natural parameterization
% ParamPost: cell (possibly empty) with natural parameters
% D: Dimenisonality of Gaussian

if (isempty(paramPost))
    nu     = randn(D,1); % Posterior mean
    Lambda = eye(D); % Posterior precision 
else % simply assigns whatever is given
    nu     = paramPost{1}; 
    Lambda = paramPost{2};
end

[ mu, Sigma, cholLambda] = getMeanFromNaturalGauss( nu, Lambda );

end



