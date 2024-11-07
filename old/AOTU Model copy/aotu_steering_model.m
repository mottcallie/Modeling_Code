% aotu_steering_model
%
% Simulates the steering behavior of a fly based on AOTU019 and AOTU025 neural input without feedforward control.
% The model generates visual object position, rotational velocity, and inputs from various neurons, 
% and computes the corresponding output based on the AOTU019 synapse type (inhibitory or excitatory).
%
% INPUTS:
%   noiTuning        - Struct containing visual tuning curves for different neurons (AOTU019, AOTU025, sum).
%                      It specifies how these neurons respond to object positions in the azimuthal space.
%   noiseLevel       - Noise level added to the rotational velocity signal.
%   startPos         - Initial object position (azimuthal angle) in degrees.
%   AOTU019synapse   - Type of synapse for AOTU019 (e.g., 'inhibitory', 'excitatory').
%   simDuration      - Duration of the simulation in seconds.
%   k                - Gain factor for computing rotational velocity from neural inputs.
%   alpha            - Alpha parameter for the activation function (ELU) of the DNa02 neuron.
%   runSettings      - Struct containing additional model parameters.
%
% OUTPUTS:
%   timebase         - Time vector for the simulation, from 0 to simDuration, based on the simulation rate (fs).
%   visobj_history   - History of the visual object positions over time for each simulation run.
%   input_history    - History of input values to the DNa02 neurons for each simulation run, for right and left sides.
%   rotvel_history   - History of the fly's rotational velocity (in degrees per second) for each simulation run.
%
% DESCRIPTION:
% This function models the behavior of a Drosophila fly during steering in response to a visual object in the absence 
% of feedforward control. It simulates the fly's interaction with the object based on inputs from visual tuning curves 
% of the AOTU019 and AOTU025 neurons, along with summed neural activity (noiSum). Depending on the type of AOTU019 synapse 
% ('inhibitory' or 'excitatory'), the model computes the input to the DNa02 motor neurons, which drive rotational velocity.
%
% The function includes delays for each neuron (AOTU019, AOTU025, and DNa02), applies visual tuning curves, and simulates 
% motor output that drives the fly's rotational velocity. Noise is added to simulate real-world variability. The rotational 
% velocity influences the perceived object position, and the function keeps track of the object's trajectory over time.
%
% The function outputs the timebase, visual object position history, inputs to DNa02 neurons, and rotational velocity 
% history, which can be used for further analysis or visualization.
%
% Example usage:
%   [timebase, visobj_history, input_history, rotvel_history] = aotu_steering_model(noiTuning, 0.05, 0, 'inhibitory', 10, 0.5, 1.0, runSettings);

function [timebase, visobj_history, input_history, rotvel_history] = aotu_steering_model(...
    noiTuning, noiseLevel, startPos, AOTU019synapse, simDuration, k, alpha, runSettings)

%% Initialize
% Set default values if any input is empty
if isempty(noiseLevel)
    noiseLevel = 0;  % Default noise level
end
if isempty(startPos)
    startPos = 0;  % Default starting position
end
if isempty(AOTU019synapse)
    AOTU019synapse = 'inhibitory';  % Default synapse type
end
if isempty(simDuration)
    simDuration = 10;  % Default simulation duration in seconds
end
if isempty(k)
    k = runSettings.k;  % Use k from runSettings if not provided
end
if isempty(alpha)
    alpha = runSettings.alpha;  % Use alpha from runSettings if not provided
end

% Extract model parameters from runSettings
visobj_position = runSettings.visObjPosition; % Object position (azimuthal space)

numRuns = runSettings.numRuns;   % Number of simulation runs
fs = runSettings.fs;             % Simulation update rate (Hz)
fpass = runSettings.fpass;       % Lowpass filter cutoff (Hz)

% Compute the number of time steps based on the simulation duration and fs
nTime = round(simDuration * fs); % Total number of time steps based on duration and update rate

%% Generate delays for each neuron in the model
delay_AOTU019 = round(runSettings.AOTU019_delay / 1000 * fs); % Delay for AOTU019 in samples
delay_Others = round(runSettings.Others_delay / 1000 * fs);   % Delay for AOTU025 and noiSum in samples
delay_DNa02 = round(runSettings.DNa02_delay / 1000 * fs);     % Delay for DNa02 output in samples

%% Generate input-output for each neuron
% Visual tuning curves for each neuron type over 360 deg azimuthal space
AOTU019R = noiTuning.AOTU019;  % Right side for AOTU019
AOTU019L = flip(AOTU019R);     % Left side (flipped)

AOTU025R = noiTuning.AOTU025;  % Right side for AOTU025
AOTU025L = flip(AOTU025R);     % Left side (flipped)

noiSumR = noiTuning.sum;       % Right side for other neurons
noiSumL = flip(noiSumR);       % Left side (flipped)

% Motor input-output for DNa02
DNa02input = runSettings.DNa02input;  % Input range for downstream neurons
DNa02output = adjELU(DNa02input, alpha, runSettings.shift); % Output range for downstream neurons

%% Initialize simulation outputs
input_history = zeros(numRuns, nTime, 2); % Inputs to downstream neurons (right/left)
DNa02R_history = zeros(numRuns, nTime);   % DNa02 right side activity
DNa02L_history = zeros(numRuns, nTime);   % DNa02 left side activity
rotvel_history = zeros(numRuns, nTime);   % Fly's rotational velocity over time

% Random noise component for rotational velocity (modeling the influence of noise)
rng(13); % Set seed for reproducibility
inputNoise = randn(numRuns, nTime);       % Random noise drawn from Gaussian distribution
inputNoise = lowpass(inputNoise', fpass, fs)'; % Low-pass filter the noise
inputNoise = noiseLevel * zscore(inputNoise);  % Scale the noise

% Initialize visual object position and history
startSide = ones(numRuns, 1); 
startSide(2:2:end) = -1;  % Alternate which side the object starts on
visobj_history = zeros(numRuns, nTime);  % Visual object position history
visobj_history(:, 1:2) = repmat((startPos * startSide), 1, 2); % Set initial position

%% Run the simulation for each time step
for t = 3:nTime
    % Compute delayed indices for AOTU019, AOTU025, and noiSum
    t_AOTU019 = max(t - delay_AOTU019, 1); % Apply delay for AOTU019
    t_Others = max(t - delay_Others, 1);   % Apply delay for AOTU025 and noiSum

    % Get previous object positions at the delayed times
    [~, p1_AOTU019] = ismember(wrapTo180(round(visobj_history(:, t_AOTU019))), visobj_position);
    [~, p2_Others] = ismember(wrapTo180(round(visobj_history(:, t_Others))), visobj_position);

    % Compute steering drive based on AOTU019 synapse type
    switch AOTU019synapse
        case "inhibitory"
            current_inputR = AOTU025R(p2_Others) - AOTU019L(p1_AOTU019) + noiSumR(p2_Others);
            current_inputL = AOTU025L(p2_Others) - AOTU019R(p1_AOTU019) + noiSumL(p2_Others);
        case "excitatory"
            current_inputR = AOTU025R(p2_Others) + AOTU019R(p1_AOTU019) + noiSumR(p2_Others);
            current_inputL = AOTU025L(p2_Others) + AOTU019L(p1_AOTU019) + noiSumL(p2_Others);
    end

    % Store the input history (right and left inputs)
    input_history(:, t, 1) = current_inputR; % Right input
    input_history(:, t, 2) = current_inputL; % Left input

    %% Apply 150 ms delay for DNa02 output
    t_DNa02 = max(t - delay_DNa02, 1);  % Apply DNa02 delay for motor output
    % Interpolate DNa02 inputs to get activity levels
    dr = interp1(DNa02input, 1:length(DNa02input), current_inputR, 'nearest');
    DNa02R_history(:, t_DNa02) = DNa02output(dr); % Right side DNa02 activity (delayed)
    dl = interp1(DNa02input, 1:length(DNa02input), current_inputL, 'nearest');
    DNa02L_history(:, t_DNa02) = DNa02output(dl); % Left side DNa02 activity (delayed)

    % Compute the rotational velocity as the difference between right and left drives
    rotvel_history(:, t) = k * (DNa02R_history(:, t_DNa02) - DNa02L_history(:, t_DNa02));

    % Update the visual object position based on rotational velocity
    visobj_history(:, t) = visobj_history(:, t) + visobj_history(:, t-1) - rotvel_history(:, t);
    
    % Add noise to the object trajectory
    visobj_history(:, t) = visobj_history(:, t) + inputNoise(:, t);
    
    % Ensure that the object position stays within the range [-180, 180]
    visobj_history(:, t) = wrapTo180(visobj_history(:, t));
end

% Generate timebase for the simulation based on actual sampling rate
timebase = linspace(0, simDuration, nTime); % Simulate from 0 to total simulation duration

end
