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

function [mdlOut] = sequential_lm(R,P1,P2,P3)
%% initialize

% fit model to 80%, test on 20%
% divide up into 5 columnns (each representing 20% of the data)
R = reshape(R,[],5);
P1  = reshape(P1,[],5);
P2  = reshape(P2,[],5);
% third predictor optional, so check first
if size(P3)<2
    nPredict = 2;
else
    nPredict = 3;
    P3  = reshape(P3,[],5);
end

% generate data storage tbl
colNames = {'Intercept','X1','R2'};
mdlOut = array2table(zeros(nPredict,3),'VariableNames',colNames);
intercept = zeros(nPredict,5);
x1 = zeros(nPredict,5);
r2 = zeros(nPredict,5);

%% fit linear model

for f = 1:5
    %figure; tiledlayout(3,1,"TileSpacing","compact")
    % set test vs fit indices
    testIdx = f;
    fitIdx = setdiff(1:5, testIdx);

    Rtrain = reshape(R(:,fitIdx),[],1);
    P1train = reshape(P1(:,fitIdx),[],1);
    P2train = reshape(P2(:,fitIdx),[],1);
    Rtest = reshape(R(:,testIdx),[],1);
    P1test = reshape(P1(:,testIdx),[],1);
    P2test = reshape(P2(:,testIdx),[],1);
    if nPredict==3
        P3train = reshape(P3(:,fitIdx),[],1);
        P3test = reshape(P3(:,testIdx),[],1);
    end

    % for each predictor
    for p = 1:nPredict
        % set model variables with each run
        switch p
            case 1
                Y = Rtrain; % response
                X = P1train; % predictor
                Xtest = P1test; % predictor for testing
            case 2
                Y = Y - Yresid; % response residual
                Rtest = Rtest-Ytest;
                X = P2train; % predictor
                Xtest = P2test; % predictor for testing
            case 3
                Y = Y - Yresid; % response residual
                Rtest = Rtest-Ytest;
                X = P3train; % predictor
                Xtest = P3test; % predictor for testing
        end

        % fit response/residual (Y) to predictor (X)
        mdl = fitlm(X,Y);

        % store model coefficients
        intercept(p,f) = mdl.Coefficients.Estimate(1);
        x1(p,f) = mdl.Coefficients.Estimate(2);
        % assess variance explained against withheld data (R2)
        Ytest = predict(mdl, Xtest);
        r2(p,f) = calculate_r_squared(Rtest,Ytest);

        % assess residual
        Yresid = predict(mdl,X);
        
        %nexttile; hold on; plot(Y); plot(X); plot(Yresid); legend({'Response';'Predictor';'Fit'})
    end
end

% output median fits
mdlOut.Intercept = median(intercept,2);
mdlOut.X1 = median(x1,2);
mdlOut.R2 = median(r2,2);

end