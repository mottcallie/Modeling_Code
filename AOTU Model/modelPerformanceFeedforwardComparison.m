% modelPerformanceFeedforwardComparison
%
% This function runs and compares the performance of a fly steering model with varying strengths of feedforward control
% using a range of steering gain (k) values. It simulates fly behavior, calculates performance metrics, and plots 
% the results to compare the effects of different feedforward scaling factors on the model's response to visual stimuli.
%
% INPUTS:
%   predicted_RF      - Struct containing the predicted receptive fields (RF) for neurons in the model.
%                       These RFs define the neural responses to visual stimuli in the azimuthal space.
%
% DESCRIPTION:
% The function initializes settings, runs the steering model across a range of feedforward scaling factors 
% (`scaleFactor`), and calculates several performance metrics (probability, variance, ISE, IAE, direction change times) 
% for each combination of steering gain values (`k`) and feedforward scaling. The `aotu_steering_model_FF` function is 
% called for each condition, where the `scaleFactor` determines the strength of the feedforward input (0 means no feedforward).
% The results are stored and plotted to visually compare how different feedforward scaling factors affect model performance.
% Summary plots for performance metrics and binned averages are also generated.
%
% STEPS:
% 1. Initialize simulation settings and set up the range of `k` values for steering gain.
% 2. Plot the receptive fields (RF) and ELU activation function for reference.
% 3. Loop over each `k` value and run the simulation for each feedforward scaling factor (from no feedforward to stronger feedforward effects).
% 4. For each run, calculate the performance metrics and store the results.
% 5. Plot example runs for selected `k` values to visually compare the behavior across different feedforward scaling factors.
% 6. Generate and save summary plots for performance metrics and binned averages.
%
% OUTPUTS:
%   The function generates and saves the following:
%   - A plot showing the predicted receptive fields and ELU function.
%   - A plot showing example runs for selected `k` values, comparing the effects of different feedforward scaling factors.
%   - Summary plots for performance metrics (probability, variance, ISE, IAE) and binned averages for object position and velocity.
%
% Example usage:
%   modelPerformanceFeedforwardComparison(predicted_RF);
%
% Note:
% The function requires additional helper functions like `aotu_steering_model_FF`, `calculatePerformanceMetrics`, 
% `analyzeObjectAndVelocityCrossings`, `analyzeErrorVsTurn`, and plotting functions such as `generateSummaryPlots` 
% and `plotRFsTogether`.

function modelPerformanceFeedforwardComparison(predicted_RF)
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
scaleFactors = 0:1:5;  % Scaling factors for the feedforward input (0 = no feedforward, 5 = max feedforward)
nScales = length(scaleFactors);  % Number of scaling factors

% Update comparison labels for feedforward scaling
comparisonLabel = arrayfun(@(x) ['ScaleFactor=' num2str(x)], scaleFactors, 'UniformOutput', false);

%% Plot RFs and ELU for reference
figure; set(gcf, 'Position', [100 100 900 400]);  % Set figure size
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile
plotRFsTogether(predicted_RF, runSettings)
nexttile
plot_adjELU(alpha, runSettings)
sgtitle('Feedforward Scaling Comparison')

% Save the plot
saveas(gcf, fullfile(folder.final, 'Feedforward_Model_Settings.png'));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, 'Feedforward_Model_Settings.svg'));

%% Run model

% Preallocate to store performance metrics
metrics_prob = zeros(nK, nScales); 
metrics_var = zeros(nK, nScales);
metrics_ISE = zeros(nK, nScales);
metrics_IAE = zeros(nK, nScales);
metrics_dirChangeTime = zeros(nK, nScales);

% Preallocate to store binned averages
obj_binned_avgs = []; 
rotvel_binned_avgs = [];
evt = [];

% Select up to X evenly spaced k values for plotting
maxPlotK = min(6, nK);  % Limit to X `k` runs for plotting
selected_k_indices = round(linspace(1, nK, maxPlotK));  % Evenly spaced indices from `kValues`

% Initialize figure for example plots
figure; set(gcf, 'Position', [100 100 1500 900]);  % Set figure size
tiledlayout(maxPlotK, nScales, 'TileSpacing', 'compact'); 
condColors = lines(nScales);  % Use a colormap to assign unique colors for each scaling factor

% Loop over k values (steering gain) for analysis
for kIdx = 1:nK
    thisK = kValues(kIdx);
    disp(['Running simulations for k = ', num2str(thisK)]);  % Display progress

    % Loop over each scaling factor (0 to 5) for the feedforward input
    for scaleIdx = 1:nScales
        scaleFactor = scaleFactors(scaleIdx);  % Current scaling factor for feedforward

        % Adjust tuning for the current scaling factor
        thisTuning = predicted_RF;

        % Run the AOTU steering model with the current scaling factor
        [timebase, visobj_history, input_history, rotvel_history] = ...
            aotu_steering_model_FF(thisTuning, noiseLevel, startPos, 'inhibitory', simDuration, thisK, alpha, scaleFactor, runSettings);

        % Calculate performance metrics
        metrics_results = calculatePerformanceMetrics(visobj_history, rotvel_history, timebase, plotSettings);
        % Store metrics for this scaling factor and k value
        metrics_prob(kIdx, scaleIdx) = metrics_results.prob;
        metrics_var(kIdx, scaleIdx) = metrics_results.var;
        metrics_ISE(kIdx, scaleIdx) = metrics_results.ISE;
        metrics_IAE(kIdx, scaleIdx) = metrics_results.IAE;

        % Measure performance at direction changes and store binned averages
        nTest = 500;
        [rotvel_cross_times, obj_bins, rotvel_bins, obj_binned_avg, rotvel_binned_avg] = analyzeObjectAndVelocityCrossings(visobj_history, rotvel_history, timebase, nTest);
        metrics_dirChangeTime(kIdx, scaleIdx) = mean(rotvel_cross_times, 'omitnan');

        % Analyze object position vs velocity
        [posvang, posBins] = analyzeErrorVsTurn(visobj_history, rotvel_history, runSettings); 

        % Store the binned averages (z represents the scaling factor dimension)
        obj_binned_avgs(kIdx, :, scaleIdx) = obj_binned_avg;
        rotvel_binned_avgs(kIdx, :, scaleIdx) = rotvel_binned_avg;
        evt(kIdx, :, scaleIdx) = posvang;

        % Plot example run
        if any(kIdx == selected_k_indices)
            % Call the function to remove large jumps
            visobj_plot = remove_large_jumps(visobj_history(1, :), 180);

            % Plot the modified visobj history
            nexttile;
            plot(timebase, visobj_plot, 'Color', condColors(scaleIdx, :));
            title(['k = ' num2str(thisK) ', ScaleFactor = ' num2str(scaleFactor)]);
            axis tight;
            if kIdx == 1
                %legend(comparisonLabel{scaleIdx});
            elseif kIdx ==nK
                xlabel('Time (s)');
            end
            if scaleIdx == 1
                ylabel('Pos')
            end
            ylim([-180 180]);
        end

    end
end

% Save the example plot comparing feedforward scaling across selected k values
sgtitle('Example Runs for Different Feedforward Scaling Factors');
saveas(gcf, fullfile(folder.final, 'Feedforward_Scaling_ExampleRuns.png'));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, 'Feedforward_Scaling_ExampleRuns.svg'));

%% Generate summary plots
% After running the model and storing all metrics and binned averages

% Define the labels for the different scaling factors
comparisonLabel = arrayfun(@(x) ['S=' num2str(x)], scaleFactors, 'UniformOutput', false);

% Call the generateSummaryPlots2 function
generateSummaryPlots2(kValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, ...
    obj_binned_avgs, rotvel_binned_avgs, evt, obj_bins, rotvel_bins, posBins, nK, comparisonLabel, folder);


end
