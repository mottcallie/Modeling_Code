% modelPerformance
%
% This function runs and compares the performance of a fly steering model across different comparison types:
% 1. Synapse type comparison (inhibitory vs excitatory)
% 2. Strength comparison (strength = 1 vs strength = 0)
% 3. Speed comparison (fast vs slow visuomotor delays)
%
% It simulates fly behavior across a range of steering gain (k) values and computes performance metrics for each 
% comparison condition. The results are stored and visualized to compare how the model behaves under different conditions.
%
% INPUTS:
%   predicted_RF      - Struct containing the predicted receptive fields (RF) for neurons in the model.
%                       These RFs define the neural responses to visual stimuli in the azimuthal space.
%   comparisonType    - Type of comparison to perform: 'synapse', 'strength', or 'speed'.
%                       - 'synapse' compares inhibitory vs excitatory synapse types.
%                       - 'strength' compares two different strength levels for the model.
%                       - 'speed' compares fast and slow visuomotor delays.
%
% DESCRIPTION:
% The function initializes settings, runs the steering model for each condition (synapse type, strength, or speed), 
% and calculates several performance metrics (probability, variance, ISE, IAE, direction change times) across a range 
% of steering gain values (`k`). The results are plotted and saved for further comparison.
%
% STEPS:
% 1. Initialize settings and determine comparison labels based on the selected `comparisonType`.
% 2. Plot the predicted receptive fields (RF) and ELU activation function for reference.
% 3. Loop over a range of `k` values (steering gain) and simulate the model for both comparison conditions.
% 4. For each run, calculate performance metrics and store the results.
% 5. Plot example runs for selected `k` values to visually compare the behavior across comparison conditions.
% 6. Generate summary plots for performance metrics and binned averages, comparing the two conditions.
%
% OUTPUTS:
%   The function generates and saves the following:
%   - A plot showing the predicted receptive fields and ELU function.
%   - A plot showing example runs for selected `k` values, comparing the two conditions.
%   - Summary plots for performance metrics (probability, variance, ISE, IAE) and binned averages for object position and velocity.
%
% Example usage:
%   modelPerformance(predicted_RF, 'synapse');  % Compare inhibitory vs excitatory synapse types
%   modelPerformance(predicted_RF, 'strength'); % Compare strength = 1 vs strength = 0
%   modelPerformance(predicted_RF, 'speed');    % Compare fast vs slow visuomotor delays
%
% Note:
% The function requires additional helper functions like `aotu_steering_model`, `calculatePerformanceMetrics`, 
% `analyzeObjectAndVelocityCrossings`, `analyzeErrorVsTurn`, and plotting functions such as `generateSummaryPlots` 
% and `plotRFsTogether`.

function modelPerformance(predicted_RF, comparisonType)
%% Initialize

% Refresh settings
[folder, plotSettings, runSettings] = modelSettings();
close all;

% Define the range to test
alpha = 0.5;
kValues = 8:2:20;  % Steering gain (k) values to iterate over
nK = length(kValues);  % Number of k values
noiseLevel = 2.5;
startPos = 0;
simDuration = 30;

% Convert alpha value to a plot-compatible string
alphaStr = sprintf('a%s', strrep(num2str(alpha, '%.2f'), '.', ''));  % E.g., 0.5 becomes a05

conditionColors = {"#0072BD"; "#D95319"};

% Variables depending on comparison type
if strcmp(comparisonType, 'synapse')
    comparisonLabel = {'Inhibitory', 'Excitatory'};
    thisSynapse = {"inhibitory"; "excitatory"};
    strengthValues = 1;  % Fixed for synapse comparison
elseif strcmp(comparisonType, 'strength')
    comparisonLabel = {'Strength 1', 'Strength 0'};
    thisSynapse = "inhibitory";  % Fixed synapse type for strength comparison
    strengthValues = [1, 0];
else  % speed comparison
    comparisonLabel = {'Fast', 'Slow'};
    thisSynapse = "inhibitory";  % Fixed synapse type for speed comparison
    visuomotorDelays = [80, 100];  % Fast = 80 ms, Slow = 100 ms
end

%% Plot RFs and ELU for reference
figure; set(gcf, 'Position', [100 100 900 400]);  % Set figure size
tiledlayout(1,2, 'TileSpacing', 'compact'); 
nexttile
plotRFsTogether(predicted_RF,runSettings)
nexttile
plot_adjELU(alpha, runSettings)
sgtitle(alphaStr)

% Save the plot
plotType = strjoin(comparisonLabel, 'v');
saveas(gcf, fullfile(folder.final, [plotType '_Model_Settings' '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, [plotType '_Model_Settings' '.svg']));

%% Run model

% Preallocate to store performance metrics and binned averages
metrics_prob = zeros(nK, 2);  % Two columns depending on comparison
metrics_var = zeros(nK, 2);
metrics_ISE = zeros(nK, 2);
metrics_IAE = zeros(nK, 2);
metrics_dirChangeTime = zeros(nK, 2);

% Preallocate to store binned averages
obj_binned_avg_1 = [];
obj_binned_avg_2 = [];
rotvel_binned_avg_1 = [];
rotvel_binned_avg_2 = [];
evt_1 = [];
evt_2 = [];

visobj_history_cond1 = [];
visobj_history_cond2 = [];

% Select up to X evenly spaced k values for plotting
maxPlotK = min(6, nK);  % Limit to X `k` runs for plotting
selected_k_indices = round(linspace(1, nK, maxPlotK));  % Evenly spaced indices from `kValues`

% Initialize figure for example plots
figure; set(gcf, 'Position', [100 100 1500 900]);  % Set figure size
tiledlayout(maxPlotK, 2, 'TileSpacing', 'compact'); 

% Loop over k values (steering gain) for analysis
for kIdx = 1:nK
    thisK = kValues(kIdx);
    % Output progress
    disp(['Running simulations for k = ', num2str(thisK)]);

    % Loop over comparison values (either synapse types, strength, or speed)
    for idx = 1:2
        if strcmp(comparisonType, 'synapse')
            thisType = thisSynapse{idx};
            thisStrength = strengthValues;  % Fixed strength for synapse comparison
        elseif strcmp(comparisonType, 'strength')
            thisStrength = strengthValues(idx);
            thisType = thisSynapse;  % Fixed synapse type for strength comparison
        else  % speed comparison
            thisStrength = 1;  % Fixed strength for speed comparison
            thisType = thisSynapse;  % Fixed synapse type
            runSettings.AOTU019_delay = visuomotorDelays(idx);
        end

        % Adjust tuning for strength if applicable
        thisTuning = predicted_RF;
        thisTuning.AOTU019 = thisTuning.AOTU019 .* thisStrength;  % Adjust tuning by strength

        % Run the AOTU steering model
        [timebase, visobj_history, input_history, rotvel_history] = ...
            aotu_steering_model(thisTuning, noiseLevel, startPos, thisType, simDuration, thisK, alpha, runSettings);

        % Calculate performance metrics
        metrics_results = calculatePerformanceMetrics(visobj_history, rotvel_history, timebase, plotSettings);
        % Store metrics for this comparison and k value
        metrics_prob(kIdx, idx) = metrics_results.prob;
        metrics_var(kIdx, idx) = metrics_results.var;
        metrics_ISE(kIdx, idx) = metrics_results.ISE;
        metrics_IAE(kIdx, idx) = metrics_results.IAE;

        % Measure performance at direction changes and store binned averages
        nTest = 500;
        [rotvel_cross_times, obj_bins, rotvel_bins, obj_binned_avgs, rotvel_binned_avgs] = analyzeObjectAndVelocityCrossings(visobj_history, rotvel_history, timebase, nTest);
        metrics_dirChangeTime(kIdx, idx) = mean(rotvel_cross_times, 'omitnan');

        % Analyze object position vs velocity
        [posvang, posBins] = analyzeErrorVsTurn(visobj_history, rotvel_history, runSettings); 

        % Store the binned averages separately for each comparison
        if idx == 1
            obj_binned_avg_1(kIdx, :) = obj_binned_avgs;
            rotvel_binned_avg_1(kIdx, :) = rotvel_binned_avgs;
            evt_1(kIdx, :) = posvang;
        else
            obj_binned_avg_2(kIdx, :) = obj_binned_avgs;
            rotvel_binned_avg_2(kIdx, :) = rotvel_binned_avgs;
            evt_2(kIdx, :) = posvang;
        end

        % Plot example run
        if any(kIdx == selected_k_indices)
            nexttile;
            % Call the function to remove large jumps
            visobj_plot = remove_large_jumps(visobj_history(1, :), 180);

            plot(timebase, visobj_plot, 'Color', conditionColors{idx});
            xlabel('Time (s)');
            ylabel('Pos');
            title(['k = ' num2str(thisK)]);
            axis tight;
            if kIdx == 1
                legend(comparisonLabel{idx});
            end
            ylim([-180 180]);

            if idx==1
                visobj_history_cond1(end+1,:) = visobj_history(1, :);
            else
                visobj_history_cond2(end+1,:) = visobj_history(1, :);
            end
        end
    end
end

% Save the example plot comparing the two conditions across selected k values
sgtitle([' Example Runs for ' comparisonLabel{1} ' vs ' comparisonLabel{2}]);
saveas(gcf, fullfile(folder.final, [plotType '_ExampleRuns' '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, [plotType '_ExampleRuns' '.svg']));

% Construct fictive path
%construct_2d_path(visobj_history_cond1(1,:), visobj_history_cond2(1,:))

%% Generate summary plots
generateSummaryPlots(kValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, ...
    obj_binned_avg_1, obj_binned_avg_2, rotvel_binned_avg_1, rotvel_binned_avg_2, ...
    evt_1, evt_2, obj_bins, rotvel_bins, posBins, nK, comparisonLabel, folder);

end
