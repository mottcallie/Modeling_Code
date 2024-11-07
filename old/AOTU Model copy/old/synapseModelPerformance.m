function synapseModelPerformance(alpha,predicted_RF)
%% Initialize

% Refresh settings
[folder, plotSettings, runSettings] = modelSettings();
close all;

% Define the range to test
kValues = 8:2:20;  % Steering gain (k) values to iterate over
nK = length(kValues);  % Number of k values
strengthValues = 1;
noiseLevel = 2.5;
startPos = 0;
simDuration = 30;
thisSynapse = {"inhibitory"; "excitatory"};

% Convert alpha value to a plot-compatible string
alphaStr = sprintf('a%s', strrep(num2str(alpha, '%.2f'), '.', ''));  % E.g., 0.5 becomes a05

%% Plot RFs and ELU for reference
figure; set(gcf, 'Position', [100 100 900 400]);  % Set figure size
tiledlayout(1,2, 'TileSpacing', 'compact'); 
nexttile
plotRFsTogether(predicted_RF,runSettings)
nexttile
plot_adjELU(alpha, runSettings)
sgtitle(alphaStr)
% Save the plot
saveas(gcf, fullfile(folder.final, ['InhibitoryvExcitatory_Model_Settings_' alphaStr '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, ['InhibitoryvExcitatory_Model_Settings_' alphaStr '.svg']));

%% Run model

% Preallocate to store performance metrics and binned averages
metrics_prob = zeros(nK, 2);  % Two columns for inhibitory and excitatory
metrics_var = zeros(nK, 2);
metrics_ISE = zeros(nK, 2);
metrics_IAE = zeros(nK, 2);
metrics_dirChangeTime = zeros(nK, 2);

% Preallocate to store binned averages
obj_binned_avg_inhibitory = [];
obj_binned_avg_excitatory = [];
rotvel_binned_avg_inhibitory = [];
rotvel_binned_avg_excitatory = [];
evt_inhibitory = [];
evt_excitatory = [];

% Select up to X evenly spaced k values for plotting
maxPlotK = min(6, nK);  % Limit to X `k` runs for plotting
selected_k_indices = round(linspace(1, nK, maxPlotK));  % Evenly spaced indices from `kValues`
% Initialize figure for example plots
figure; set(gcf, 'Position', [100 100 1500 900]);  % Set figure size
tiledlayout(maxPlotK, 2, 'TileSpacing', 'compact');
synapseColors = {"#0072BD";"#D95319"};

% Loop over k values (steering gain) for analysis
for kIdx = 1:nK
    thisK = kValues(kIdx);
    % Output progress
    disp(['Running simulations for k = ', num2str(thisK)]);

    % Loop over synapse types (inhibitory and excitatory)
    for synapseIdx = 1:2
        thisSynapseType = thisSynapse{synapseIdx};

        % Run the AOTU steering model with adjusted synapse type
        [timebase, visobj_history, input_history, rotvel_history] = ...
            aotu_steering_model(predicted_RF, noiseLevel, startPos, thisSynapseType, simDuration, thisK, alpha, runSettings);

        % Calculate performance metrics
        metrics_results = calculatePerformanceMetrics(visobj_history, rotvel_history, timebase, plotSettings);
        % Store metrics for this synapse type and k value
        metrics_prob(kIdx, synapseIdx) = metrics_results.prob;
        metrics_var(kIdx, synapseIdx) = metrics_results.var;
        metrics_ISE(kIdx, synapseIdx) = metrics_results.ISE;
        metrics_IAE(kIdx, synapseIdx) = metrics_results.IAE;

        % Measure performance at direction changes and store binned averages
        nTest = 100;
        [rotvel_cross_times, obj_bins, rotvel_bins, obj_binned_avgs, rotvel_binned_avgs] = analyzeObjectAndVelocityCrossings(visobj_history, rotvel_history, timebase, nTest);
        % Store direction changes
        metrics_dirChangeTime(kIdx, synapseIdx) = mean(rotvel_cross_times, 'omitnan');

        % Analyze object position v velocity
        [posvang, posBins] = analyzeErrorVsTurn(visobj_history, rotvel_history, runSettings); 

        % Store the binned averages for inhibitory and excitatory separately
        if strcmp(thisSynapseType, 'inhibitory')
            obj_binned_avg_inhibitory(kIdx, :) = obj_binned_avgs;
            rotvel_binned_avg_inhibitory(kIdx, :) = rotvel_binned_avgs;
            evt_inhibitory(kIdx, :) = posvang;
        else
            obj_binned_avg_excitatory(kIdx, :) = obj_binned_avgs;
            rotvel_binned_avg_excitatory(kIdx, :) = rotvel_binned_avgs;
            evt_excitatory(kIdx, :) = posvang;
        end

        % Plot example run
        if any(kIdx==selected_k_indices)
            nexttile;
            plot(timebase, visobj_history(1, :),'Color',synapseColors{synapseIdx});
            xlabel('Time (s)');
            ylabel('Pos');
            title(['k = ' num2str(thisK)]);
            axis tight;
            if kIdx==1
                legend(thisSynapseType)
            end
            ylim([-180 180]);
        end
    end
end

% Save the example plot
sgtitle([(alphaStr) ' Example Runs for Inhibitory vs Excitatory']);
saveas(gcf, fullfile(folder.final, ['InhibitoryvExcitatory_ExampleRuns_' alphaStr '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, ['InhibitoryvExcitatory_ExampleRuns_' alphaStr '.svg']));

%% Generate summary plots
% Summary Plot Comparing Across Different Metrics (Inhibitory vs Excitatory)
figure;
set(gcf, 'Position', [100 100 400 800]);  % Set figure size
tiledlayout(4, 1, 'TileSpacing', 'compact');  % One row for each metric

% Plot 1: Probability near setpoint across k values
nexttile;
plot(kValues, metrics_prob(:, 1), '-o', 'DisplayName', 'Inhibitory');
hold on;
plot(kValues, metrics_prob(:, 2), '-o', 'DisplayName', 'Excitatory');
xlabel('Steering Gain (k)');
ylabel('Probability Near Setpoint');
legend('show');
title('Probability Near Setpoint (Inhibitory vs Excitatory)');
axis padded
ylim([0 0.5])
grid on;

% Plot 2: Circular variance across k values
nexttile;
plot(kValues, metrics_var(:, 1), '-o', 'DisplayName', 'Inhibitory');
hold on;
plot(kValues, metrics_var(:, 2), '-o', 'DisplayName', 'Excitatory');
xlabel('Steering Gain (k)');
ylabel('1- Circular Variance');
title('Circular Variance (Inhibitory vs Excitatory)');
axis padded
ylim([0 1])
grid on;

% Plot 3: Integral of Squared Error (ISE) across k values
nexttile;
plot(kValues, metrics_ISE(:, 1), '-o', 'DisplayName', 'Inhibitory');
hold on;
plot(kValues, metrics_ISE(:, 2), '-o', 'DisplayName', 'Excitatory');
xlabel('Steering Gain (k)');
ylabel('ISE');
title('Squared Error (Inhibitory vs Excitatory)');
axis padded
grid on;

% Plot 4: Integral of Absolute Error (IAE) across k values
nexttile;
plot(kValues, metrics_IAE(:, 1), '-o', 'DisplayName', 'Inhibitory');
hold on;
plot(kValues, metrics_IAE(:, 2), '-o', 'DisplayName', 'Excitatory');
xlabel('Steering Gain (k)');
ylabel('IAE');
title('Absolute Error (Inhibitory vs Excitatory)');
axis padded
grid on;

% Save the summary plot comparing inhibitory and excitatory across metrics
sgtitle({[alphaStr 'Summary of Model Performance'],'Across Metrics (Inhibitory vs Excitatory)'});
saveas(gcf, fullfile(folder.final, ['InhibitoryvExcitatory_Compare_SummaryMetrics_' alphaStr '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, ['InhibitoryvExcitatory_Compare_SummaryMetrics_' alphaStr  '.svg']));

%% Generate two summary plots for binned averages vs. object position and rotational velocity
% Call the generateSummaryPlots function after the model run loop
comparisonLabel = {'Inhibitory', 'Excitatory'};

% Generate the summary plots
generateSummaryPlots(kValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, ...
    obj_binned_avg_inhibitory, obj_binned_avg_excitatory, rotvel_binned_avg_inhibitory, rotvel_binned_avg_excitatory, ...
    evt_inhibitory, evt_excitatory, obj_bins, rotvel_bins, posBins, nK, alphaStr, comparisonLabel, folder);

end