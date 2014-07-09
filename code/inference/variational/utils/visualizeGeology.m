function  h = visualizeGeology(geologySettings, sensorGravLocations, theta2Transitions, ctrlPts, y )
%VISUALIZEGEOLOGY Summary of this function goes here
%   y: Observations
transitionFunction = @(query)theta2Transitions(query, ctrlPts);
h = figure;
visGravWorld(geologySettings,transitionFunction,sensorGravLocations, y);
subplot(1,2,1); title('Density');  colorbar; view(3);
subplot(1,2,2);  title('Gravity Observations'); colorbar;
 

end

