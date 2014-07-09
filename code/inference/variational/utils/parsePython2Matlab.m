function  str  = parsePython2Matlab( str )
%PARSEPYTHON2MATLAB Summary of this function goes here
%   Detailed explanation goes here

str = strrep(str, 'np.', '');
str = strrep(str, 'f.shape', 'size(theta)');
str = strrep(str, 'f', 'theta');
str = strrep(str, '**', '.^');


return;



