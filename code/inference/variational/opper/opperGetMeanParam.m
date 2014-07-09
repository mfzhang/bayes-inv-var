function [ mq, Sigmaq Lambdaq] = opperGetMeanParam( omegaq, eta, Sigmap, invSigmap )
%OPPERGETMEANPARAM Get mean parameters from Opper parameterization
%   Detailed explanation goes here
% mq     = (Sigmap)(omegaq) --> omegaq = (Sigmap{^-1})(mq)
% Sigmaq = (diag(eta) + So^{-1})^{-1} 

mq            = Sigmap*omegaq;
Lambdaq       = diag(eta) + invSigmap; % precision
cholLambdaq   = getChol(Lambdaq);  
Sigmaq        = getInverseChol(cholLambdaq);


return;

