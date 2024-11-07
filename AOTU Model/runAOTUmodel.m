% runAOTUmodel
%
% Main script for running multiple iterations of a simplified AOTU (anterior optic tubercle) model 
% to evaluate the impact of varying parameters such as synapse weight, synapse sign (inhibitory vs excitatory),
% noise levels, open-loop vs closed-loop conditions, and feedforward modulation on steering behavior.
%
% Created: 10/05/2024 - MC
% Updated: 11/01/2024 - MC
%
%% Initialize workspace and variables
% Clear all variables and close any open figures to ensure a clean workspace
clear
close all

% Initialize folder paths, plot settings, and run-specific parameters by calling the modelSettings function.
[folder, plotSettings, runSettings] = modelSettings();

% Change the directory to the script's file path to access data files.
cd(folder.filePath)

% Load the receptive field (RF) data from the 'Pursuit_RFs.mat' file, which contains information about
% the RFs of AOTU neurons for visual pursuit behavior.
load("Pursuit_RFs.mat");

% Process the raw RF data using the 'processRFdata' function to generate predicted receptive fields (predicted_RF).
% These predicted RFs will be used in subsequent model performance comparisons.
predicted_RF = processRFdata(Pursuit_RFs);

%% Calibrate alpha parameter
modelCalibration(predicted_RF)

% Plot RFs and ELU for reference
plotModelSettings(predicted_RF, runSettings)

% Save the plot
saveas(gcf, fullfile(folder.settings, ['Model_Settings' '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, ['Model_Settings' '.svg']));

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
% Compare modrel performance with and without feedforward modulation from the AOTU019 neuron.
% This comparison evaluates the influence of feedforward input on the fly's ability to maintain stable pursuit behavior.
modelPerformance(predicted_RF, 'feedforward')
close all  % Close all open figures after the comparison is complete.
