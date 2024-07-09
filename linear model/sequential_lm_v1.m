% function for performing a sequential linear model
% Briefly, the response data will be fit to the first predictor. The
% residual will then be fit to the second predictor, and so forth. Variance
% explained will be assessed by calculating the R2 value at each step.
%
% INPUT
% R - Response array
% P1 - First predictor array
% P2 - Second predictor array
% P3 - Third predictor array (optional)
%
% OUTPUT
% mdlOut - Variance explained
%
% CREATED 04/22/24 MC
%

function [mdlOut] = sequential_lm_v1(R,P1,P2,P3)
%% initialize

% check number of trials in total
nTrial = size(R,2);

% check how many predictors are being tested
if size(P3)<2
    nPredict = 2;
else
    nPredict = 3;
end

% generate data storage tbl
colNames = {'Intercept','X1','R2'};
mdlOut = array2table(zeros(nPredict,3),'VariableNames',colNames);

%% fit linear model

% for each predictor
for p = 1:nPredict
    % set model variables with each run
    switch p
        case 1
            Y = R; % response
            X = P1; % first predictor
        case 2
            Y = Y - yFit; % response residual
            X = P2; % first predictor
        case 3
            Y = Y - yFit; % response residual
            X = P3; % first predictor
    end

    % fit response/residual to predictor
    mdl = fitlm(X,Y);

    % store model coefficients
    mdlOut.Intercept(p) = mdl.Coefficients.Estimate(1);
    mdlOut.X1(p) = mdl.Coefficients.Estimate(2);
    % store variance explained (R2)
    mdlOut.R2(p) = mdl.Rsquared.Adjusted;

    % generate model prediction
    yFit = predict(mdl,X);

end

end