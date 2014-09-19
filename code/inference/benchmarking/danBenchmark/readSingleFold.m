%% Read data
function [xtrain, ftrain, ytrain, xtest, ftest, ytest] = ...
               readSingleFold(x, f, y, train, test, k)
    xtrain = x(train(k,:),:);
    ftrain = f(train(k,:),:);
    ytrain = y(train(k,:),:);

    xtest = x(test(k,:),:);
    ftest = f(test(k,:),:);
    ytest = y(test(k,:),:);
    
return;


