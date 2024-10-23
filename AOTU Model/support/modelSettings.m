% modelSettings
%
% This function initializes and returns the folder paths, plot settings, and run settings used in the steering model simulations.
% It sets up essential parameters required for running the fly steering simulations, including simulation parameters,
% neuron-related settings, and paths for saving output data and plots.
%
% OUTPUTS:
%   folder        - A struct containing paths to directories for saving various types of data.
%   plotSettings  - A struct containing settings for plotting the model outputs.
%   runSettings   - A struct containing model-specific parameters for the simulation.
%
% DESCRIPTION:
% The `modelSettings` function initializes key parameters and settings that are used across different simulation runs and visualizations.
% It defines parameters for the visual object position, neuron delays, low-pass filter settings for noise, and the receptive field colors.
% Additionally, the function defines folder paths where simulation outputs and visualizations will be saved, and configures general plotting properties.
%
% Example usage:
%   [folder, plotSettings, runSettings] = modelSettings();

function [folder, plotSettings, runSettings] = modelSettings()
%% Run settings (model-specific parameters)
% Visual object positions in azimuthal space
runSettings.visObjPosition = linspace(-180, 180, 361);

% Simulation parameters
runSettings.numRuns = 20000; % Number of times to run the simulation (ideal > 20k)
runSettings.k = 20;          % Max change in head direction per time step
runSettings.fs = 100;        % Simulation update rate, in Hz
runSettings.fpass = 2;       % Lowpass cutoff for random component, in Hz

% Neuron settings
runSettings.AOTU019_delay = 80; %ms
runSettings.Others_delay = 100; %ms
runSettings.DNa02_delay = 150; %ms

% Open-loop parameters
runSettings.waveAmp = 35;    % Amplitude of the wave, in degrees

% Downstream process parameters
runSettings.minIn = -1;      % Minimum input value for downstream process
runSettings.maxIn = 3.5;     % Maximum input value for downstream process
runSettings.alpha = 0.075;     % alpha, specifies slope of negative non-linear portion
runSettings.shift = 0.025;    % shift, determines at what point linear/non-linear portions begin
runSettings.DNa02input = linspace(runSettings.minIn, runSettings.maxIn, 1000); % Input range for DNa02
runSettings.DNa02output = adjELU(runSettings.DNa02input, runSettings.alpha, runSettings.shift); % Output curve based on ELU nonlinearity

%% Folder paths
folder.thisFile = matlab.desktop.editor.getActiveFilename;
[folder.filePath,~,~] = fileparts(folder.thisFile); % Get file path of current script
folder.trials = [folder.filePath '/trials'];        % Path for storing trial data
folder.vectors = [folder.filePath '/vectors'];        % Path for storing vectors
folder.summary = [folder.filePath '/summary'];      % Path for storing summary data
folder.final = [folder.filePath '/final'];      % Path for storing final data

%% Plot settings
plotSettings.rfTicks = -180:60:180; % Angle ticks for receptive field plotting (from -180 to 180 degrees)
plotSettings.nEx = 5;               % Number of example model traces to plot

% Plotting - velocity parameters
plotSettings.vMax = 75;             % Maximum velocity for binning
plotSettings.vBinSize = 5;          % Size of velocity bins
plotSettings.vEdges = -plotSettings.vMax+plotSettings.vBinSize/2:plotSettings.vBinSize:plotSettings.vMax; % Velocity bin edges
plotSettings.vCentr = plotSettings.vEdges(1:end-1)+plotSettings.vBinSize/2; % Centers of velocity bins

% Plotting - heading and position parameters
plotSettings.nPolarBins = 20;       % Number of bins for polar histograms (used for heading distribution)
plotSettings.hdBinSize = 10;        % Bin size for standard histograms (used for position data)

% General plot settings
plotSettings.dim = [100 100 1800 800]; % Window dimensions for the plot
plotSettings.celltypesMain = {}; % Placeholder for RF cell types to be populated later
plotSettings.rfColors = {"#77AC30";"#0072BD";"#7E2F8E"}; % Colors for RF curves (sum, AOTU019, AOTU025)
plotSettings.vHistRange = [0 0.2];  % Range for velocity histograms
plotSettings.hHistRange = [0 0.4];  % Range for heading histograms
plotSettings.opacity = 0.4;         % Opacity level for plot shading

end