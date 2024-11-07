% analyzeErrorVsTurn
% Determine the relationship between setpoint error (object position) and the fly's turn behavior.
% Bins velocity by object position, combining all trials into a single array.
%
% INPUTS:
%   visobj_history  - Object position history (numRuns x nTime).
%   rotvel_history  - Rotational velocity history (numRuns x nTime).
%
% OUTPUTS:
%   posvang         - Binned angular velocity vs object position (averaged across trials).
%   posBins         - Object position bins.
%
function [posvang, posBins] = analyzeErrorVsTurn(visobj_history, rotvel_history,runSettings)
    %% Initialize
    % Set bin parameters
    posMax = 180;  % +/- deg
    posBin = 2;   % deg

    % Create bins
    posEdge = -posMax - posBin / 2 : posBin : posMax + posBin / 2;  % Bin edges for object position
    posBins = -posMax : posBin : posMax;                            % Bin centers for object position

    % Initialize array for binned angular velocities
    posvang = nan(length(posBins), 1);
    minBin = 500;  % Minimum number of points required per bin

    %% (optional) shift according to lag estimates
    % fetch shift indices for each lag
    lag_ms = mean([runSettings.AOTU019_delay,runSettings.Others_delay]);
    lag_idx = round(lag_ms / 1000 * runSettings.fs);

    % shift and exclude data at start/stop of trial
    rotvel_history = circshift(rotvel_history,-lag_idx,2); %shift
    rotvel_history(:,end-lag_idx+1:end) = nan; %exclude shifts

    %% Bin angular velocity by object position
    % Combine all trials into a single array
    % Reshape and combine all trials into a single array
    thisVisObj = reshape(visobj_history, [], 1);
    thisRotVel = reshape(rotvel_history, [], 1);

    % Discretize object position
    discPanelps = discretize(thisVisObj, posEdge, posBins);

    for p = 1:length(posBins)
        thisBinIdx = find(discPanelps == posBins(p));
        if length(thisBinIdx) >= minBin
            posvang(p) = mean(thisRotVel(thisBinIdx), 'omitnan');
        end
    end
end