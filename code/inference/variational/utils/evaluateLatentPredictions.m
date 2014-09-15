function evaluateLatentPredictions(xtrain, ftrain, xtest, muStar, sigmaStar, ftest, strTitle)
%% Evaluate Latent Predicitons
figure;
if size(xtest,1)
    [foo idx] = sort(xtest);
    plot_confidence_interval(xtest(idx),muStar(idx),sigmaStar(idx),1.96);
    hold on;
    plot(xtest(idx), ftest(idx), 'g--', 'MarkerSize', 12);
    plot(xtrain, ftrain, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    xlabel('x');
    ylabel('f');
    title(['Latent Function: ', strTitle ], 'FontSize', 18);
    legend({'Pred Error Bars', 'Pred Mean', 'True Test', 'Training'});
    
else
    disp('This problem is not 1-dimesional, plots not shown');
end

return;

