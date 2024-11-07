% aotu_steering_model_FF_OL
%
% Simulates the steering behavior of a fly using a feedforward model for the AOTU019 neuron in an open-loop configuration. 
% The model incorporates visual object position, rotational velocity, and inputs from various neurons, and generates corresponding
% outputs based on feedforward control scaled by a given factor. The visual object follows a predefined triangle wave trajectory.
%
% INPUTS:
%   noiTuning        - Struct containing visual tuning curves for various neurons (AOTU019, AOTU025, sum).
%                      It defines the response of these neurons to object positions in the azimuthal space.
%   noiseLevel       - Noise level to add to the rotational velocity signal.
%   startPos         - Starting object position (azimuthal angle) in degrees.
%   AOTU019synapse   - Type of synapse for AOTU019 (e.g., 'inhibitory', 'excitatory').
%   simDuration      - Duration of the simulation in seconds.
%   k                - Gain factor for computing rotational velocity from neural inputs.
%   alpha            - Alpha parameter used in the activation function (ELU) for the DNa02 neuron output.
%   scaleFactor      - Scaling factor applied to the feedforward term (0 means no feedforward).
%   runSettings      - Struct containing additional simulation parameters (e.g., number of runs, sample rate, delays, etc.).
%
% OUTPUTS:
%   timebase         - Time vector for the simulation, from 0 to simDuration, based on the sampling rate (fs).
%   visobj_history   - History of the visual object positions over time for each simulation run. The object follows 
%                      a triangle wave trajectory in open-loop configuration.
%   input_history    - History of input values to the downstream DNa02 neurons for each simulation run, for right and left sides.
%   rotvel_history   - History of the fly's rotational velocity (in degrees per second) for each simulation run.
%
% DESCRIPTION:
% This function models the behavior of a Drosophila fly during steering in response to a visual object in an open-loop configuration.
% The visual object follows a predefined alternating triangle wave pattern, simulating oscillatory motion. The fly's rotational velocity 
% is computed based on inputs from visual tuning curves of different neurons (AOTU019, AOTU025) and summed neuron activity (noiSum).
% Feedforward control is incorporated based on the fly's rotational velocity and scaled by the input `scaleFactor`. When `scaleFactor`
% is 0, there is no feedforward effect.
%
% The model incorporates delays for each neuron type (AOTU019, AOTU025, DNa02), applies visual tuning curves to compute inputs to the 
% DNa02 neurons, and simulates the corresponding motor output, which drives the fly's rotational velocity. Noise is added to the signal 
% to simulate real-world variability in the fly's behavior. The fly's perceived object position is influenced by its rotational velocity, 
% and noise is applied to the visual object trajectory as well.
% 
% The function outputs the timebase, visual object position history (in open-loop), inputs to the DNa02 neurons, and rotational velocity 
% history for further analysis or plotting.
%
% Example usage:
%   [timebase, visobj_history, input_history, rotvel_history] = aotu_steering_model_FF_OL(noiTuning, 0.05, 0, 'inhibitory', 10, 0.5, 1.2, runSettings);
%
function [timebase, visobj_history, input_history, rotvel_history] = aotu_steering_model_FF(...
    noiTuning, noiseLevel, startPos, AOTU019synapse, simDuration, k, alpha, scaleFactor, runSettings)

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
    t_DNa02 = max(t - delay_DNa02, 1);  % Apply DNa02 delay for motor output

    % Get previous object positions at the delayed times
    [~, p1_AOTU019] = ismember(wrapTo180(round(visobj_history(:, t_AOTU019))), visobj_position);
    [~, p2_Others] = ismember(wrapTo180(round(visobj_history(:, t_Others))), visobj_position);

    % Compute steering drive based on AOTU019 synapse type
    % Feedforward effect (scaling by scaleFactor)
    turning_right = rotvel_history(:, t_AOTU019) * scaleFactor;
    turning_right(turning_right<0) = 0;  % Only positive right turns contribute
    turning_left = -rotvel_history(:, t_AOTU019) * scaleFactor;
    turning_left(turning_left<0) = 0;    % Only positive left turns contribute

    % Compute input values based on feedforward scaling
    current_inputR = AOTU025R(p2_Others) - (AOTU019L(p1_AOTU019) + turning_left) + noiSumR(p2_Others);
    current_inputL = AOTU025L(p2_Others) - (AOTU019R(p1_AOTU019) + turning_right) + noiSumL(p2_Others);
    
    % Ensure inputs do not exceed limits
    current_inputR(current_inputR < -1) = -1;
    current_inputL(current_inputL < -1) = -1;

    % Store the input history (right and left inputs)
    input_history(:, t, 1) = current_inputR; % Right input
    input_history(:, t, 2) = current_inputL; % Left input

    %% Apply delay for DNa02 output
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
