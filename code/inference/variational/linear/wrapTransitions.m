function theta2Transitions  = wrapTransitions( query, ctrlPts, ...
                                                    transitions, fixedParams, dimCtrl )
%WRAPTRANSITIONS Summary of this function goes here
%   Detailed explanation goes here

C = reshape(ctrlPts,dimCtrl,dimCtrl);
theta2Transitions =  transitions(query, [fixedParams C ]);

return;





