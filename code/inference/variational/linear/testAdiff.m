% tests adiff object
clear all; clc;


addpath(genpath('~/Dropbox/Matlab/autodiff'));


x       = rand(5,1);
A       = rand(10,5);
funcPtr = @(xx)  A*xx;

adObj   = adiff(x);
fwdObj  = feval(funcPtr, adObj);

val = adiffget(fwdObj,'value')
J = adiffget(fwdObj,'derivative')
