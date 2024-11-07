% modelCalibration
%
% This function runs and compares the performance of a fly steering model
% across a range of steering gain (k) values and sensitivity levels (alpha).
% It simulates fly behavior under each combination of k and alpha, computes
% performance metrics, and visualizes the results.
%
% INPUTS:
%   predicted_RF - Struct containing the predicted receptive fields (RF) for
%                  neurons in the model. These RFs define the neural responses
%                  to visual stimuli in the azimuthal space.
%
% DESCRIPTION:
% The function initializes model settings, runs the steering model across 
% specified k and alpha values, and calculates several performance metrics, 
% including:
%   - Circular Variance: Measures consistency of the model's position.
%   - Probability Near Setpoint: Probability of staying within +/- 5 degrees.
%   - Integral of Squared Error (ISE): Measures squared error relative to setpoint.
%   - Integral of Absolute Error (IAE): Measures absolute error relative to setpoint.
%
% For each k and alpha combination, the function stores and visualizes metrics:
%   1. 1 - Variance plotted across gains (k) for each alpha.
%   2. Input-output curve showing downstream neuron response (DNa02 output) as a
%      function of the input and alpha.
%
% The results are saved as PNG and SVG files in designated folders for further analysis.
%
% CREATED: 11/01/2024 - MC
%
function modelCalibration(predicted_RF)
%% Initialize

% Refresh settings
[folder, plotSettings, runSettings] = modelSettings();
close all;

% Define the range to test for alpha and k
runSettings.numRuns = 100; %override
alpha = [0.05, 0.075, 0.1, 0.25, 0.5, 1];  % Range of alpha values
nA = length(alpha);  % Number of alpha values
kValues = 4:2:18;  % Range of k values (steering gain)
nK = length(kValues);  % Number of k values
noiseLevel = 1.5;
startPos = 0;
simDuration = 30;

%% Run model
% Preallocate to store performance metrics
metrics_var = zeros(nK, nA);

% Loop over k values (steering gain) for analysis
for kIdx = 1:nK
    thisK = kValues(kIdx);
    disp(['Running simulations for k = ', num2str(thisK)]);

    % Loop over alpha values
    for alphaIdx = 1:nA
        thisAlpha = alpha(alphaIdx);

        % Output progress for alpha
        disp(['  Alpha = ', num2str(thisAlpha)]);

        % Run the AOTU steering model for this combination of k and alpha
        [timebase, visobj_history, ~, rotvel_history] = aotu_steering_model(predicted_RF, noiseLevel, startPos, "inhibitory", simDuration, thisK, thisAlpha, runSettings);

        % Calculate performance metrics
        metrics_results = calculatePerformanceMetrics(visobj_history, rotvel_history, timebase, plotSettings);

        % Store metrics for this combination of k and alpha
        metrics_var(kIdx, alphaIdx) = metrics_results.var;
        
        % Additional metrics can also be stored similarly, e.g., metrics_prob, metrics_ISE, etc.
    end
end


% Plot results
% Create a tiled layout with two rows
figure
tiledlayout(1,2);

% First Tile: 1 - Variance across gains for each alpha
nexttile;
hold on;

% Generate the x-axis values from kValues
x = kValues;

% Plot 1 - variance for each alpha condition
for alphaIdx = 1:nA
    % Calculate 1 - variance for each k value for the current alpha
    y = 1 - metrics_var(:, alphaIdx);
    
    % Plot the line for this alpha, using a different color for each line
    plot(x, y, 'DisplayName', ['\alpha = ' num2str(alpha(alphaIdx))]);
end

% Customize the plot for 1 - variance
xlabel('Steering Gain (k)');
ylabel('1 - Variance');
title('Model Performance');
legend('show');  % Displays a legend for each alpha value
grid on;
hold off;

% Second Tile: Input-Output Curve for each Alpha
nexttile;
hold on;

% Loop over each alpha to calculate and plot input-output curves
for alphaIdx = 1:nA
    % Define input and calculate output using adjELU
    DNa02input = runSettings.DNa02input;  % Input range for downstream neurons
    DNa02output = adjELU(DNa02input, alpha(alphaIdx), runSettings.shift); % Output for each alpha
    
    % Plot input vs. output for this alpha value
    plot(DNa02input, DNa02output, 'DisplayName', ['\alpha = ' num2str(alpha(alphaIdx))]);
end

% Customize the plot for input-output curve
xlabel('DNa02 Input');
ylabel('DNa02 Output');
title('ELU');
legend('show');  % Displays a legend for each alpha value
grid on;
hold off;

% Save the example plot comparing the two conditions across selected k values
saveas(gcf, fullfile(folder.settings, ['Model_Calibration' '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, ['Model_Calibration' '.svg']));

end
