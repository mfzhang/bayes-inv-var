
function visGravWorld(geologySettings,transitionFunction,sensorLocations,sensorMeasurements)


qGrid = wbt.prior.query3D([geologySettings.boundaries(1)...
    geologySettings.boundaries(3) geologySettings.boundaries(5)],...
    [geologySettings.boundaries(2),geologySettings.boundaries(4),...
    geologySettings.boundaries(6)], [100 100]);
gridVals = transitionFunction(qGrid);


subplot(1,2,1);
regionColours = {'r','g','b',0.5*[1 1 1],'y'};
regionAlpha = [.15 .15 .15 .8 .15];
wbt.vis.plotWorldGeometryVol(geologySettings, qGrid, gridVals,regionColours, regionAlpha);
xlabel('longitude'); ylabel('latitude'); zlabel('depth')
view(0,90)
drawnow


subplot(1,2,2);
scatter3(sensorLocations(1,:),sensorLocations(2,:),...
    sensorMeasurements,20,sensorMeasurements,'filled')
view(0,90);title('Sensor Readings')
xlabel('longitude'); ylabel('latitude'); zlabel('Gravity Anomaly')
axis square
