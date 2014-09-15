function evaluateObservablePredictions(xtrain, xtest, ytest, yStar, strTitle)
if (size(xtrain,2) == 1)
    [foo idx] = sort(xtest);
    xtest = xtest(idx);
    ytest = ytest(idx);
    yStar = yStar(idx);
    figure;
    plot(xtest, yStar, 'k-', 'LineWidth', 2, 'MarkerSize', 12);
    hold on;
    plot(xtest, ytest, 'g--', 'MarkerSize', 12);
    hold on;
    plot(xtrain, min(ytest)*ones(size(xtrain,1),1), 'rx', 'MarkerSize', 12);

    title(['Observable Function:', strTitle ], 'FontSize', 18);
    legend({'Predicted', 'True', 'Training Input Location'});
end

return;



