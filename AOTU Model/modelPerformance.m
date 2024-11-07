% modelPerformance
%
% This function runs and compares the performance of a fly steering model
% across different comparison types. It simulates fly behavior across a
% range of steering gain (k) values and computes performance metrics for each 
% comparison condition. The results are stored and visualized to compare how
% the model behaves under different conditions.
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
% CREATED: 10/30/2024 - MC
%

function modelPerformance(predicted_RF, comparisonType)
%% Initialize

% Refresh settings
[folder, plotSettings, runSettings] = modelSettings();
close all;

% Define the range to test
alpha = 0.5;
kValues = 2:2:12;  % Steering gain (k) values to iterate over
nK = length(kValues);  % Number of k values
noiseLevel = 1;
startPos = 0;
simDuration = 30;
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
elseif strcmp(comparisonType, 'speed')
    comparisonLabel = {'Fast', 'Slow'};
    thisSynapse = "inhibitory";  % Fixed synapse type for speed comparison
    visuomotorDelays = [runSettings.AOTU019_delay, runSettings.Others_delay];  % Fast = 80 ms, Slow = 100 ms
elseif strcmp(comparisonType, 'feedforward')
    comparisonLabel = {'repeat019', 'feedforward'};
    thisSynapse = {"feedforwardctrl"; "feedforward"};
end

% Run model
% Preallocate to store performance metrics and binned averages
metrics_prob = zeros(nK, 2);  % Two columns depending on comparison
metrics_var = zeros(nK, 2);
metrics_ISE = zeros(nK, 2);
metrics_IAE = zeros(nK, 2);
metrics_dirChangeTime = zeros(nK, 2);

% Preallocate to store binned averages
indist_1 = [];
indist_2 = [];
rotvel_binned_avg_1 = [];
rotvel_binned_avg_2 = [];
objvel_binned_avg_1 = [];
objvel_binned_avg_2 = [];
evt_1 = [];
evt_2 = [];

visobj_history_cond1 = [];
visobj_history_cond2 = [];

steering_drive = [];

% Select up to X evenly spaced k values for plotting
maxPlotK = min(6, nK);  % Limit to X `k` runs for plotting
selected_k_indices = round(linspace(1, nK, maxPlotK));  % Evenly spaced indices from `kValues`

% Initialize figure for example plots
figure; set(gcf, 'Position', [100 100 1500 900]);  % Set figure size
tiledlayout(maxPlotK, 2, 'TileSpacing', 'compact'); 

% Loop over k values (steering gain) for analysis
for kIdx = 1:nK
    thisK = kValues(kIdx);
    runSettings.nK = kIdx;
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
        elseif  strcmp(comparisonType, 'speed')
            thisStrength = 1;  % Fixed strength for speed comparison
            thisType = thisSynapse;  % Fixed synapse type
            runSettings.AOTU019_delay = visuomotorDelays(idx);
        elseif strcmp(comparisonType, 'feedforward')
            thisStrength = 1;
            thisType = thisSynapse{idx};
        end

        % Adjust tuning for strength if applicable
        thisTuning = predicted_RF;
        thisTuning.AOTU019 = thisTuning.AOTU019 .* thisStrength;  % Adjust tuning by strength

        % Run the AOTU steering model
        [timebase, visobj_history, input_history, rotvel_history] = aotu_steering_model(thisTuning, noiseLevel, startPos, thisType, simDuration, thisK, alpha, runSettings);

        % Calculate performance metrics
        metrics_results = calculatePerformanceMetrics(visobj_history, rotvel_history, timebase, plotSettings);
        % Store metrics for this comparison and k value
        metrics_prob(kIdx, idx) = metrics_results.prob;
        metrics_var(kIdx, idx) = metrics_results.var;
        metrics_ISE(kIdx, idx) = metrics_results.ISE;
        metrics_IAE(kIdx, idx) = metrics_results.IAE;

        % Measure performance at direction changes and store binned averages
        nTest = runSettings.numRuns;
        [rotvel_cross_times, rotvel_bins, rotvel_binned_avgs,objvel_binned_avgs] = analyzeObjectAndVelocityCrossings(visobj_history, rotvel_history, timebase, nTest);
        metrics_dirChangeTime(kIdx, idx) = mean(rotvel_cross_times, 'omitnan');

        % Analyze object position vs velocity
        [posvang, posBins] = analyzeErrorVsTurn(visobj_history, rotvel_history, runSettings); 
        
        % Compute normalized distribution of input history
        [input_distribution, inBins] = compute_normalized_distribution(input_history, runSettings);
        % Store the binned averages separately for each comparison
        if idx == 1
            rotvel_binned_avg_1(kIdx, :) = rotvel_binned_avgs;
            objvel_binned_avg_1(kIdx, :) = objvel_binned_avgs;
            evt_1(kIdx, :) = posvang;
            indist_1(kIdx,:) = input_distribution;
        else
            rotvel_binned_avg_2(kIdx, :) = rotvel_binned_avgs;
            objvel_binned_avg_2(kIdx, :) = objvel_binned_avgs;
            evt_2(kIdx, :) = posvang;
            indist_2(kIdx,:) = input_distribution;
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
        
        if kIdx==1
            % Check relationship between input and steering
            %steering_drive(idx,:) = aotu_steering_model_position(thisTuning, thisType, alpha, runSettings);
        end
    end
end

% Save the example plot comparing the two conditions across selected k values
sgtitle([' Example Runs for ' comparisonLabel{1} ' vs ' comparisonLabel{2}]);
saveas(gcf, fullfile(folder.final, [comparisonLabel{1} 'v' comparisonLabel{2} '_ExampleRuns' '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, [comparisonLabel{1} 'v' comparisonLabel{2} '_ExampleRuns' '.svg']));

%% Generate summary plots
% Generate and save summary metrics plot
generateSummaryMetrics(kValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, comparisonLabel, folder);

% Generate and save summary plot of binned crossing times
generateSummaryCrossing(nK, kValues, rotvel_binned_avg_1, rotvel_binned_avg_2, objvel_binned_avg_1, objvel_binned_avg_2, rotvel_bins, comparisonLabel, folder);

% Generate and save summary plot of Angular Velocity vs Object Position (EVT)
generateSummaryEVT(nK, kValues, evt_1, evt_2, posBins, comparisonLabel,folder);

% Generate and save summary plot for steering drive
%generateSteeringPlot(steering_drive, comparisonLabel, runSettings,folder)

% Generate and save summary plot for inputs to DNa02
generateSummaryInputDistributions(nK, kValues, indist_1, indist_2, inBins, comparisonLabel, runSettings, folder)

end
