% runAOTUmodel
%
% Main script for running multiple iterations of a simplified AOTU (anterior optic tubercle) model 
% to evaluate the impact of varying parameters such as synapse weight, synapse sign (inhibitory vs excitatory),
% noise levels, open-loop vs closed-loop conditions, and feedforward modulation on steering behavior.
%
% The script runs several comparison models using the AOTU019 neuron with different conditions, including:
% - Alpha values affecting the activation function.
% - Synapse strength comparison (active vs silenced).
% - Synapse sign comparison (inhibitory vs excitatory).
% - Visuomotor delay comparison (fast vs normal).
% - Feedforward input comparison (with vs without feedforward).
%
% Each comparison evaluates the model’s performance across these conditions using preloaded receptive field (RF) data for pursuit behavior.
% Performance metrics are generated and visualized to analyze how the model reacts under different parameter settings.
%
% Created: 10/05/2024 - MC
%
%% Initialize workspace and variables
% Clear all variables and close any open figures to ensure a clean workspace
clear
close all

% Initialize folder paths, plot settings, and run-specific parameters by calling the modelSettings function.
% This sets up essential configurations needed for the simulations.
[folder, plotSettings, runSettings] = modelSettings();

% Change the directory to the script's file path to access data files.
cd(folder.filePath)

% Load the receptive field (RF) data from the 'Pursuit_RFs.mat' file, which contains information about
% the RFs of AOTU neurons for visual pursuit behavior.
load("Pursuit_RFs.mat");

% Process the raw RF data using the 'processRFdata' function to generate predicted receptive fields (predicted_RF).
% These predicted RFs will be used in subsequent model performance comparisons.
predicted_RF = processRFdata(Pursuit_RFs);

%% Compare model performance across different alpha values
% This section runs model performance comparisons across a range of alpha values, which modify the activation function.
% Alpha values control the non-linear portion of the DNa02 neuron output, affecting the model's response to visual input.

% Define the alpha values to test in the model simulations.
alpha = [0.1, 0.25, 0.5, 0.8, 1]; 
nAlpha = length(alpha);  % Calculate the number of alpha values to iterate over.

% Loop over each alpha value and run the model performance comparison.
for aIdx = 1:nAlpha
    thisAlpha = alpha(aIdx);  % Current alpha value for the iteration
    disp(['   Preparing batch for alpha = ' num2str(thisAlpha)])  % Display progress in the command window.
    
    % Call 'compareModelPerformance' to evaluate model performance for the current alpha value.
    % This function generates performance metrics based on the steering behavior of the model with the specified alpha.
    compareModelPerformance(thisAlpha, predicted_RF)
end

% Close all open figures to avoid memory issues and overlapping visualizations.
close all

%% Run AOTU019 active vs silenced comparison model
% Compare model performance when AOTU019 is fully active vs silenced (strength = 1 vs strength = 0).
% This comparison evaluates how the presence or absence of AOTU019 neuron activity influences steering behavior.
modelPerformance(predicted_RF, 'strength');  % Specify 'strength' as the comparison type.
close all  % Close all open figures after the comparison is complete.

%% Run AOTU019 inhibitory vs excitatory comparison model
% Compare model performance when AOTU019 uses inhibitory synapses vs excitatory synapses.
% This comparison evaluates how the sign of the synapse (inhibitory or excitatory) affects the fly’s steering response.
modelPerformance(predicted_RF, 'synapse');  % Specify 'synapse' as the comparison type.
close all  % Close all open figures after the comparison is complete.

%% Run AOTU019 fast vs normal comparison model
% Compare model performance between fast and normal visuomotor delay conditions (80 ms vs 100 ms).
% This comparison assesses how changes in the visuomotor delay affect steering accuracy and timing.
modelPerformance(predicted_RF, 'speed');  % Specify 'speed' as the comparison type.
close all  % Close all open figures after the comparison is complete.

%% Run AOTU019 with vs without feedforward input comparison model
% Compare model performance with and without feedforward modulation from the AOTU019 neuron.
% This comparison evaluates the influence of feedforward input on the fly's ability to maintain stable pursuit behavior.
modelPerformanceFeedforwardComparison(predicted_RF);  % Call the feedforward comparison function.
close all  % Close all open figures after the comparison is complete.

%% Run AOTU019 with vs without feedforward input comparison model open loop
% Compare model performance with and without feedforward modulation from the AOTU019 neuron.
% This comparison evaluates the influence of feedforward input on the fly's ability to maintain stable pursuit behavior.
modelPerformanceFeedforwardComparison2(predicted_RF);  % Call the feedforward comparison function.
close all  % Close all open figures after the comparison is complete.
